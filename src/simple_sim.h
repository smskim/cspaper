#ifndef CONFIDENCESEQUENCES_SIMPLE_SIM_H_
#define CONFIDENCESEQUENCES_SIMPLE_SIM_H_

#include <algorithm>
#include <limits>
#include <memory>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/zeta.hpp>

#include "uniform_boundaries.h"

class RunningStats {
 public:

  int t;
  double sum_X;
  double sum_X_sq;
  double sum_sq_prediction_errors;

  RunningStats() : t(0), sum_X(0), sum_X_sq(0), sum_sq_prediction_errors(1e-6)
      {}

  void add(double x) {
    double x_hat = mean();
    t += 1;
    sum_X += x;
    sum_X_sq += pow(x, 2);
    sum_sq_prediction_errors += pow(x - x_hat, 2);
  }

  double mean() const {
    return t == 0 ? 0 : sum_X / t;
  }

  double quadratic_variation(double mean) const {
    // max with epsilon to deal with numerical error
    return std::max(1e-6, sum_X_sq + t * pow(mean, 2) - 2 * mean * sum_X);
  }

  RunningStats negate() const {
    RunningStats negated_stats = *this;
    negated_stats.sum_X *= -1;
    return negated_stats;
  }

  double martingale(const double mean) const {
    return sum_X - t * mean;
  }
};


class Support {
  public:

  double min;
  double max;

  Support() : min(-std::numeric_limits<double>::infinity()),
              max(std::numeric_limits<double>::infinity()) {}
  Support(double min, double max) : min(min), max(max) {}

  Support negate() const {return Support(-max, -min);}
};

class Strategy {
  public:

  virtual ~Strategy() {}
  virtual bool can_reject(double mean, const RunningStats& stats,
                          const Support support) const = 0;
  virtual double find_lower_confidence_bound(const RunningStats& stats,
                                             const Support support) const = 0;
};

class NormalMixtureStrategy : public Strategy {
 public:

  NormalMixtureStrategy(double t_opt, double alpha)
      : t_opt_(t_opt), alpha_(alpha) {}

  bool can_reject(double mean, const RunningStats& stats,
                  const Support support) const;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const;

  static double v_factor(const Support support) {
    return pow(support.max - support.min, 2) / 4.0;
  }

 private:
  confseq::TwoSidedNormalMixture get_mixture(const Support support) const {
    return confseq::TwoSidedNormalMixture(v_factor(support) * t_opt_, alpha_);
  }

  const int t_opt_;
  const double alpha_;
};

class BetaBinomialStrategy : public Strategy {
 public:

  BetaBinomialStrategy(double v_opt, double alpha)
      : v_opt_(v_opt), alpha_(alpha) {}

  bool can_reject(double mean, const RunningStats& stats,
                  const Support support) const;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const;

 private:
  static double get_v(const int t, const double mean, const Support support) {
    return pow(support.max - support.min, 2) / 4.0;
  }

  confseq::BetaBinomialMixture get_mixture(const double mean,
                                           const Support support) const {
    assert(support.min < mean & mean < support.max);
    const double g = mean - support.min;
    const double h = support.max - mean;
    return confseq::BetaBinomialMixture(v_opt_, alpha_, g, h, false);
  }

  const int v_opt_;
  const double alpha_;
};

class EmpiricalBernsteinStrategy : public Strategy {
 public:

  EmpiricalBernsteinStrategy(double v_opt, double alpha)
      : v_opt_(v_opt), alpha_(alpha) {}

  bool can_reject(double mean, const RunningStats& stats,
                  const Support support) const;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const;

 private:
  confseq::GammaExponentialMixture get_mixture(const Support support) const {
    return confseq::GammaExponentialMixture(v_opt_, alpha_ / 2,
                                            support.max - support.min);
  }

  const int v_opt_;
  const double alpha_;
};

class NaiveSelfNormalizedStrategy : public Strategy {
 public:

  NaiveSelfNormalizedStrategy(double v_opt, double alpha)
      : mixture_(v_opt, alpha), alpha_(alpha) {}

  bool can_reject(double mean, const RunningStats& stats,
                  const Support support) const;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const;

 private:
  const confseq::TwoSidedNormalMixture mixture_;
  const double alpha_;
};

class ZTestStrategy : public Strategy {
 public:
  ZTestStrategy(const double alpha) : z_factor_(get_z_factor(alpha)) {}

  static double get_z_factor(const double alpha) {
    boost::math::normal dist(0.0, 1.0);
    return quantile(dist, 1 - alpha / 2.0);
  }

  bool can_reject(double mean, const RunningStats& stats, const Support support)
      const override;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const override;

  double z_test_radius(const int t, const double mean, const Support support)
      const;

 private:
  const double z_factor_;
};

class PointwiseBernoulliStrategy : public Strategy {
 public:
  PointwiseBernoulliStrategy(const double alpha) : alpha_(alpha) {}

  bool can_reject(double mean, const RunningStats& stats, const Support support)
      const override;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const override;

 private:
  double bernoulli_psi_star(double u, double g, double h) const;
  double tail_bound_exponent(const double mean, const RunningStats& stats,
                             const Support support) const;
  const double alpha_;
};

class PointwiseHoeffdingStrategy : public Strategy {
 public:
  PointwiseHoeffdingStrategy(const double alpha) : alpha_(alpha) {}

  bool can_reject(double mean, const RunningStats& stats, const Support support)
      const override;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const override;

 private:
  double radius(const int t, const Support support) const;
  const double alpha_;
};

class LinearBoundStrategy : public Strategy {
 public:
  LinearBoundStrategy(const int t_opt, const double alpha)
      : t_opt_(t_opt), alpha_(alpha) {}

  bool can_reject(double mean, const RunningStats& stats, const Support support)
      const override;
  double find_lower_confidence_bound(const RunningStats& stats,
                                     const Support support) const override;

 private:
  double boundary(const int t, const Support support) const;
  const int t_opt_;
  const double alpha_;
};

#endif // CONFIDENCESEQUENCES_SIMPLE_SIM_H_
