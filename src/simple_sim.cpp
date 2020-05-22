#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <future>
#include <iostream>
#include <limits>
#include <random>
#include <thread>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/roots.hpp>

#include <RcppCommon.h>

#include "simple_sim.h"

class SimulationConfig;
RCPP_EXPOSED_CLASS(SimulationConfig)

#include <Rcpp.h>

using namespace Rcpp;

using StepGenerator =
    std::function<std::function<double(std::default_random_engine&)>()>;

class DiscreteStep{
 public:
  DiscreteStep(std::vector<double> values, std::vector<double> probs)
      : values(values), distribution(probs.begin(), probs.end()) {}

  double operator()(std::default_random_engine& rng) {
    return values[distribution(rng)];
  }

 private:
  const std::vector<double> values;
  std::discrete_distribution<> distribution;
};

bool NormalMixtureStrategy::can_reject(double mean, const RunningStats& stats,
                                       const Support support) const {
  const double s = stats.martingale(mean);
  const double v = v_factor(support) * stats.t;
  return get_mixture(support).log_superMG(s, v) >= log(1 / alpha_);
}

double NormalMixtureStrategy::find_lower_confidence_bound(
    const RunningStats& stats, const Support support) const {
  const double v = v_factor(support) * stats.t;
  const double radius =
      get_mixture(support).bound(v, log(1 / alpha_)) / stats.t;
  return std::max(support.min, stats.mean() - radius);
}

bool BetaBinomialStrategy::can_reject(double mean, const RunningStats& stats,
                                      const Support support) const {
  const double s = stats.martingale(mean);
  const double v = (support.max - mean) * (mean - support.min) * stats.t;
  return get_mixture(mean, support).log_superMG(s, v) >= log(1 / alpha_);
}

double bernoulli_lower_bound(
    const double empirical_mean, const int num_trials, const Support support,
    const double alpha, const double v_opt, const double alpha_opt) {
  const double threshold = log(1 / alpha);
  const double empirical_p =
      (empirical_mean - support.min) / (support.max - support.min);
  if (empirical_p <= 1e-8) {
    return 0;
  }

  auto objective = [empirical_mean, num_trials, support, v_opt, alpha_opt,
                    threshold](const double p) {
    if (p <= 1e-8) {
      return 1.0;
    } else if (p >= 1 - 1e-8) {
      return -1.0;
    } else {
      const double mean = support.min + p * (support.max - support.min);
      const double g = mean - support.min;
      const double h = support.max - mean;
      const confseq::BetaBinomialMixture mixture(v_opt, alpha_opt, g, h, false);
      const double log_superMG = mixture.log_superMG(
          (empirical_mean - mean) * num_trials, g * h * num_trials);
      return log_superMG - threshold;
    }
  };

  assert(objective(empirical_p) < 0);
  boost::math::tools::eps_tolerance<double> tolerance(40);
  auto lower_bound_pair = boost::math::tools::bisect(
      objective, 0.0, empirical_p, tolerance);
  const double p_lower_bound =
      (lower_bound_pair.first + lower_bound_pair.second) / 2.0;
  return support.min + p_lower_bound * (support.max - support.min);
}

double BetaBinomialStrategy::find_lower_confidence_bound(
    const RunningStats& stats, const Support support) const {
  return bernoulli_lower_bound(
      stats.mean(), stats.t, support, alpha_, v_opt_, alpha_);
}

bool EmpiricalBernsteinStrategy::can_reject(
    double mean, const RunningStats& stats, const Support support) const {
  const double s = abs(stats.martingale(mean));
  const double v = stats.sum_sq_prediction_errors;
  return get_mixture(support).log_superMG(s, v) >= log(2 / alpha_);
}

double EmpiricalBernsteinStrategy::find_lower_confidence_bound(
    const RunningStats& stats, const Support support) const {
  const double radius = get_mixture(support).bound(
      stats.sum_sq_prediction_errors, log(2 / alpha_)) / stats.t;
  return std::max(support.min, stats.mean() - radius);
}

bool NaiveSelfNormalizedStrategy::can_reject(
    double mean, const RunningStats& stats, const Support support) const {
  const double s = stats.martingale(mean);
  const double v = stats.sum_sq_prediction_errors;
  return mixture_.log_superMG(s, v) >= log(1 / alpha_);
}

double NaiveSelfNormalizedStrategy::find_lower_confidence_bound(
    const RunningStats& stats, const Support support) const {
  const double radius = mixture_.bound(
      stats.sum_sq_prediction_errors, log(1 / alpha_)) / stats.t;
  return std::max(support.min, stats.mean() - radius);
}

double ZTestStrategy::z_test_radius(const int t, const double mean,
                                    const Support support) const {
  double sample_p = (mean - support.min) / (support.max - support.min);
  double std_dev = sqrt(t * sample_p * (1 - sample_p))
      * (support.max - support.min);
  return z_factor_ * std_dev;
}

bool ZTestStrategy::can_reject(double mean, const RunningStats& stats,
                               const Support support) const {
  double radius = z_test_radius(stats.t, mean, support);
  return abs(stats.martingale(mean)) >= radius;
}

double ZTestStrategy::find_lower_confidence_bound(const RunningStats& stats,
                                                  const Support support) const {
  double ci_radius = z_test_radius(stats.t, stats.mean(), support)
      / stats.t;
  return std::max(stats.mean() - ci_radius, support.min);
}

double PointwiseBernoulliStrategy::bernoulli_psi_star(double u, double g,
                                                      double h) const {
  double q = g * (1 + h * u) / (g + h);
  assert(-1e-8 <= q && q <= 1 + 1e-8);
  double p = g / (g + h);
  if (q <= 0) {
    return -log(1 - p) / (g * h);
  } else if (q >= 1) {
    return -log(p) / (g * h);
  } else {
    return 1 / (g * h) * (q * log(q / p) + (1 - q) * log((1 - q) / (1 - p)));
  }
}

double PointwiseBernoulliStrategy::tail_bound_exponent(
    const double mean, const RunningStats& stats, const Support support) const {
  const double g = mean - support.min;
  const double h = support.max - mean;
  const double s = stats.martingale(mean);
  const double v = g * h * stats.t;
  return v * bernoulli_psi_star(s / v, g, h);
}

bool PointwiseBernoulliStrategy::can_reject(
    double mean, const RunningStats& stats, const Support support) const {
  return tail_bound_exponent(mean, stats, support) >= log(2 / alpha_);
}

double PointwiseBernoulliStrategy::find_lower_confidence_bound(
    const RunningStats& stats, const Support support) const {
  const double threshold = log(2 / alpha_);
  boost::math::tools::eps_tolerance<double> tolerance(40);
  auto objective = [this, stats, support, threshold](const double mean) {
    if (mean <= support.min) {
      return 1.0;
    } else if (mean >= support.max) {
      return -1.0;
    } else {
      return tail_bound_exponent(mean, stats, support) - threshold;
    }
  };

  if (stats.mean() > support.min) {
    auto lower_bound_pair = boost::math::tools::bisect(
        objective, support.min, stats.mean(), tolerance);
    return (lower_bound_pair.first + lower_bound_pair.second) / 2.0;
  } else {
    return support.min;
  }
}

double PointwiseHoeffdingStrategy::radius(const int t, const Support support)
    const {
  return sqrt(2 * NormalMixtureStrategy::v_factor(support) * t *
              log(2 / alpha_));
}

bool PointwiseHoeffdingStrategy::can_reject(
    double mean, const RunningStats& stats, const Support support) const {
  return abs(stats.martingale(mean)) >= radius(stats.t, support);
}

double PointwiseHoeffdingStrategy::find_lower_confidence_bound(
    const RunningStats& stats, const Support support) const {
  return std::max(stats.mean() - radius(stats.t, support) / stats.t,
                  support.min);
}

double LinearBoundStrategy::boundary(const int t, const Support support) const {
  const double v_factor = NormalMixtureStrategy::v_factor(support);
  const double m = v_factor * t_opt_;
  const double Vt = v_factor * t;
  const double x = sqrt(2 * m * log(2 / alpha_));
  return x + x / (2 * m) * (Vt - m);
}

bool LinearBoundStrategy::can_reject(
    double mean, const RunningStats& stats, const Support support) const {
  return abs(stats.martingale(mean)) >= boundary(stats.t, support);
}

double LinearBoundStrategy::find_lower_confidence_bound(
    const RunningStats& stats, const Support support) const {
  return std::max(stats.mean() - boundary(stats.t, support) / stats.t,
                  support.min);
}

double get_width(const Strategy& strategy, const RunningStats& stats,
                 const Support support) {
  double lower = strategy.find_lower_confidence_bound(stats, support);
  double upper = -strategy.find_lower_confidence_bound(stats.negate(),
                                                       support.negate());
  return upper - lower;
}

class SimulationConfig {
 public:
  SimulationConfig(int num_walks, int num_steps, int num_threads,
                   std::vector<double> evaluate_width_times)
      : num_threads(num_threads), num_walks(num_walks), num_steps(num_steps),
        evaluate_width_times(evaluate_width_times) {
    std::vector<double> values = {-1, 1};
    std::vector<double> probs = {0.5, 0.5};
    discrete_steps(values, probs);
    normal_mixture(10, 0.05);
  }

  SimulationConfig& discrete_steps(std::vector<double> values,
                                   std::vector<double> probs) {
    draw_step = [=]() {return DiscreteStep(values, probs);};
    support = Support(*std::min_element(values.begin(), values.end()),
                      *std::max_element(values.begin(), values.end()));
    return *this;
  }

  SimulationConfig& normal_mixture(double t_opt, double alpha) {
    strategy = std::make_shared<NormalMixtureStrategy>(t_opt, alpha);
    return *this;
  }

  SimulationConfig& beta_binomial_mixture(double v_opt, double alpha) {
    strategy = std::make_shared<BetaBinomialStrategy>(v_opt, alpha);
    return *this;
  }

  SimulationConfig& empirical_bernstein(double v_opt, double alpha) {
    strategy = std::make_shared<EmpiricalBernsteinStrategy>(v_opt, alpha);
    return *this;
  }

  SimulationConfig& naive_self_normalized(double v_opt, double alpha) {
    strategy = std::make_shared<NaiveSelfNormalizedStrategy>(v_opt, alpha);
    return *this;
  }

  SimulationConfig& pointwise_bernoulli(double alpha) {
    strategy = std::make_shared<PointwiseBernoulliStrategy>(alpha);
    return *this;
  }

  SimulationConfig& z_test(double alpha) {
    strategy = std::make_shared<ZTestStrategy>(alpha);
    return *this;
  }

  SimulationConfig& pointwise_hoeffding(double alpha) {
    strategy = std::make_shared<PointwiseHoeffdingStrategy>(alpha);
    return *this;
  }

  SimulationConfig& linear_bound(int t_opt, double alpha) {
    strategy = std::make_shared<LinearBoundStrategy>(t_opt, alpha);
    return *this;
  }

  const int num_threads;
  const int num_walks;
  const int num_steps;
  const std::vector<double> evaluate_width_times;
  StepGenerator draw_step;
  std::shared_ptr<const Strategy> strategy;
  Support support;
};

int simulate_one_walk(const SimulationConfig& config,
                      std::default_random_engine &rng,
                      std::vector<double>& walk_widths) {
  RunningStats stats;
  int stopping_time = -1;
  auto step_fn = config.draw_step();
  auto next_width_time = config.evaluate_width_times.begin();
  for (int t = 1; t <= config.num_steps; t++) {
    double step = step_fn(rng);
    stats.add(step);

    if (*next_width_time == t) {
      walk_widths.push_back(get_width(*config.strategy, stats, config.support));
      ++next_width_time;
    }

    if (stopping_time == -1
        && config.strategy->can_reject(0, stats, config.support)) {
      stopping_time = t;
    }
  }
  return stopping_time;
}

void simulate_many_walks(const int thread_id, const SimulationConfig* config,
                         const bool* const stop_all_threads,
                         std::vector<int>* stopping_times,
                         std::vector<std::vector<double>>* widths) {
  std::default_random_engine rng;
  const int walks_per_thread = config->num_walks / config->num_threads;
  const int start_index = thread_id * walks_per_thread;

  for (int walk_index = start_index;
       walk_index < start_index + walks_per_thread;
       walk_index++) {
    rng.seed(293887823 + walk_index);
    if (*stop_all_threads) break;

    (*widths)[walk_index].reserve(config->evaluate_width_times.size());
    int stopping_time = simulate_one_walk(*config, rng, (*widths)[walk_index]);
    (*stopping_times)[walk_index] += stopping_time;

    if (thread_id == 0 && ((walk_index + 1) % 10) == 0) {
      fprintf(stderr, "%d\r", (walk_index + 1) * config->num_threads);
      fflush(stderr);
    }
  }
}

List run_simulation(const SimulationConfig& config) {
  const int actual_num_walks =
      (config.num_walks / config.num_threads) * config.num_threads;
  bool stop_all_threads = false;
  std::future<void> futures[config.num_threads];
  std::vector<int> stopping_times(actual_num_walks);
  std::vector<std::vector<double>> widths(actual_num_walks);
  auto start = std::chrono::steady_clock::now();

  for (int thread_id = 0; thread_id < config.num_threads; thread_id++) {
    futures[thread_id] = std::async(
        std::launch::async, simulate_many_walks,
        thread_id, &config, &stop_all_threads, &stopping_times, &widths);
  }

  int num_finished = 0;
  while (!stop_all_threads && num_finished < config.num_threads) {
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    try {
      Rcpp::checkUserInterrupt();
    } catch (Rcpp::internal::InterruptedException) {
      fprintf(stderr, "Caught interrupt, stopping threads\n");
      stop_all_threads = true;
      for (auto& future : futures) {
        future.wait();
      }
      throw;
    }
    num_finished = 0;
    for (auto& future : futures) {
      auto status = future.wait_for(std::chrono::milliseconds(0));
      if (status == std::future_status::ready) ++num_finished;
    }
  }
  for (auto& future : futures) future.get(); // check for exceptions

  int num_crossed = std::count_if(
      stopping_times.begin(),
      stopping_times.end(),
      [](int i){return i > 0;});
  auto end = std::chrono::steady_clock::now();
  int elapsed_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
      .count();
  fprintf(stderr, "%d / %d walks crossed, time elapsed=%.3fs\n",
          num_crossed, actual_num_walks, elapsed_ms / 1000.0);

  NumericMatrix width_matrix(
      actual_num_walks,
      config.evaluate_width_times.size()
                             );
  for (int walk_index = 0; walk_index < actual_num_walks; walk_index++) {
    std::vector<double> walk_widths = widths[walk_index];
    for (int step = 0; step < config.evaluate_width_times.size(); step++) {
      width_matrix(walk_index, step) =
          walk_widths.size() > step ? walk_widths[step] : NA_REAL;
    }
  }
  return List::create(
      _["num_walks"] = actual_num_walks,
      _["num_steps"] = config.num_steps,
      _["stopping_times"] = stopping_times,
      _["evaluate_width_times"] = config.evaluate_width_times,
      _["widths"] = width_matrix);
}

RCPP_MODULE(simple_sim_cpp) {
  class_<SimulationConfig>("SimulationConfig")
      .constructor<int, int, int, std::vector<double>>()
      .method("discrete_steps", &SimulationConfig::discrete_steps)
      .method("normal_mixture", &SimulationConfig::normal_mixture)
      .method("beta_binomial_mixture", &SimulationConfig::beta_binomial_mixture)
      .method("empirical_bernstein", &SimulationConfig::empirical_bernstein)
      .method("naive_self_normalized", &SimulationConfig::naive_self_normalized)
      .method("pointwise_bernoulli", &SimulationConfig::pointwise_bernoulli)
      .method("z_test", &SimulationConfig::z_test)
      .method("pointwise_hoeffding", &SimulationConfig::pointwise_hoeffding)
      .method("linear_bound", &SimulationConfig::linear_bound);

  function("run_simulation", &run_simulation);
}
