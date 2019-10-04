/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2019 Folkert Nobels    (nobels@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

#include <fenv.h>

/* Local headers. */
#include "swift.h"

/**
 * @brief Compute the Pearson correlation coefficient for two sets of numbers
 *
 * The pearson correlation coefficient between two sets of numbers can be
 * calculated as:
 *
 *           <x*y> - <x>*<y>
 * r_xy = ----------------------
 *         (var(x) * var(y))^.5
 *
 * In the case that both sets are purely uncorrelated the value of the
 * Pearson correlation function is expected to be close to 0. In the case that
 * there is positive correlation r_xy > 0 and in the case of negative
 * correlation, the function has r_xy < 0.
 *
 * @param mean1 average of first series of numbers
 * @param mean2 average of second series of numbers
 * @param total12 sum of x_i * y_i of both series of numbers
 * @param var1 variance of the first series of numbers
 * @param var2 variance of the second series of numbers
 * @param counter number of elements in both series
 * @return the Pearson correlation coefficient
 * */
double pearsonfunc(double mean1, double mean2, double total12, double var1,
                   double var2, int counter) {

  const double mean12 = total12 / (double)counter;
  const double correlation = (mean12 - mean1 * mean2) / sqrt(var1 * var2);
  return fabs(correlation);
}

/**
 * @brief Test to check that the pseodo-random numbers in SWIFT are random
 * enough for our purpose.
 *
 * The test initializes with the current time and than creates 20 ID numbers
 * it runs the test using these 20 ID numbers. Using these 20 ID numbers it
 * Checks 4 different things:
 * 1. The mean and variance are correct for random numbers generated by this
 *    ID number.
 * 2. The random numbers from this ID number do not cause correlation in time.
 *    Correlation is checked using the Pearson correlation coefficient which
 *    should be sufficiently close to zero.
 * 3. A small offset in ID number of 2, doesn't cause correlation between
 *    the two sets of random numbers (again with the Pearson correlation
 *    coefficient) and the mean and variance of this set is
 *    also correct.
 * 4. Different physical processes in random.h are also uncorrelated and
 *    produce the correct mean and variance as expected. Again the correlation
 *    is calculated using the Pearson correlation coefficient.
 *
 * More information about the Pearson correlation coefficient can be found in
 * the function pearsonfunc above this function.
 *
 * @param argc Unused
 * @param argv Unused
 * @return 0 if everything is fine, 1 if random numbers are not random enough.
 */
int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  uint64_t cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* Time-step size */
  const int time_bin = 30;

  /* Try a few different values for the ID */
  for (int i = 0; i < 10; ++i) {

    const int64_t id = rand() * (1LL << 31) + rand();
    const integertime_t increment = (1LL << time_bin);
    const int64_t idoffset = id + 2;

    message("Testing id=%lld time_bin=%d", id, time_bin);

    double total = 0., total2 = 0.;
    int count = 0;

    /* Pearson correlation variables for different times */
    double sum_previous_current = 0.;
    double previous = 0.;

    /* Pearson correlation for two different IDs */
    double pearsonIDs = 0.;
    double totalID = 0.;
    double total2ID = 0.;

    /* Pearson correlation for different processes */
    double pearson_star_sf = 0.;
    double pearson_star_se = 0.;
    double pearson_star_bh = 0.;
    double pearson_sf_se = 0.;
    double pearson_sf_bh = 0.;
    double pearson_se_bh = 0.;

    /* Calculate the mean and <x^2> for these processes */
    double total_sf = 0.;
    double total_se = 0.;
    double total_bh = 0.;

    double total2_sf = 0.;
    double total2_se = 0.;
    double total2_bh = 0.;

    /* Check that the numbers are uniform over the full-range of useful
     * time-steps */
    for (integertime_t ti_current = 0LL; ti_current < max_nr_timesteps;
         ti_current += increment) {

      ti_current += increment;

      const double r =
          random_unit_interval(id, ti_current, random_number_star_formation);
      if (r < 0.0 || r >= 1.0) {
        error("Generated random vlaue %f is not in [0, 1).", r);
      }

      total += r;
      total2 += r * r;
      count++;

      /* Calculate for correlation between time.
       * For this we use the pearson correlation of time i and i-1 */
      sum_previous_current += r * previous;
      previous = r;

      /* Calculate if there is a correlation between different ids */
      const double r_2ndid = random_unit_interval(idoffset, ti_current,
                                                  random_number_star_formation);
      if (r_2ndid < 0.0 || r_2ndid >= 1.0) {
        error("Generated random vlaue %f is not in [0, 1).", r_2ndid);
      }

      /* Pearson correlation for small different IDs */
      pearsonIDs += r * r_2ndid;
      totalID += r_2ndid;
      total2ID += r_2ndid * r_2ndid;

      /* Calculate random numbers for the different processes and check
       * that they are uncorrelated */

      const double r_sf =
          random_unit_interval(id, ti_current, random_number_stellar_feedback);
      if (r_sf < 0.0 || r_sf >= 1.0) {
        error("Generated random vlaue %f is not in [0, 1).", r_sf);
      }

      const double r_se = random_unit_interval(
          id, ti_current, random_number_stellar_enrichment);
      if (r_se < 0.0 || r_se >= 1.0) {
        error("Generated random vlaue %f is not in [0, 1).", r_se);
      }

      const double r_bh =
          random_unit_interval(id, ti_current, random_number_BH_feedback);
      if (r_bh < 0.0 || r_bh >= 1.0) {
        error("Generated random vlaue %f is not in [0, 1).", r_bh);
      }

      /* Calculate the correlation between the different processes */
      total_sf += r_sf;
      total_se += r_se;
      total_bh += r_bh;

      total2_sf += r_sf * r_sf;
      total2_se += r_se * r_se;
      total2_bh += r_bh * r_bh;

      pearson_star_sf += r * r_sf;
      pearson_star_se += r * r_se;
      pearson_star_bh += r * r_bh;
      pearson_sf_se += r_sf * r_se;
      pearson_sf_bh += r_sf * r_bh;
      pearson_se_bh += r_se * r_bh;
    }

    const double mean = total / (double)count;
    const double var = total2 / (double)count - mean * mean;

    /* Pearson correlation calculation for different times */
    // const double mean_xy = sum_previous_current / ((double)count - 1.f);
    // const double correlation = (mean_xy - mean * mean) / var;
    const double correlation =
        pearsonfunc(mean, mean, sum_previous_current, var, var, count - 1);

    /* Mean for different IDs */
    const double meanID = totalID / (double)count;
    const double varID = total2ID / (double)count - meanID * meanID;

    /* Pearson correlation between different IDs*/
    const double correlationID =
        pearsonfunc(mean, meanID, pearsonIDs, var, varID, count);

    /* Mean and <x^2> for different processes */
    const double mean_sf = total_sf / (double)count;
    const double mean_se = total_se / (double)count;
    const double mean_bh = total_bh / (double)count;

    const double var_sf = total2_sf / (double)count - mean_sf * mean_sf;
    const double var_se = total2_se / (double)count - mean_se * mean_se;
    const double var_bh = total2_bh / (double)count - mean_bh * mean_bh;

    /* Correlation between different processes */
    const double corr_star_sf =
        pearsonfunc(mean, mean_sf, pearson_star_sf, var, var_sf, count);
    const double corr_star_se =
        pearsonfunc(mean, mean_se, pearson_star_se, var, var_se, count);
    const double corr_star_bh =
        pearsonfunc(mean, mean_bh, pearson_star_bh, var, var_bh, count);
    const double corr_sf_se =
        pearsonfunc(mean_sf, mean_se, pearson_sf_se, var_sf, var_se, count);
    const double corr_sf_bh =
        pearsonfunc(mean_sf, mean_bh, pearson_sf_bh, var_sf, var_bh, count);
    const double corr_se_bh =
        pearsonfunc(mean_se, mean_bh, pearson_se_bh, var_se, var_bh, count);

    /* Verify that the mean and variance match the expected values for a uniform
     * distribution */

    /* Set the allowed standard deviation */
    const double std_check = 5.;

    /* The mean is expected to deviate a maximum of std_check * std / sqrt(N) */
    const double tolmean = std_check / sqrtf(12.f * count);

    /* the variance is expected to deviate a maximum of std_check * variance
     * * sqrt(2/(n-1)) */
    const double tolvar =
        std_check * sqrtf(2.f / (12.f * ((double)count - 1.f)));

    /* The correlation coefficient is expected to deviate sqrt(1-R^2)
     * / sqrt(n-2), in our case <R> = 0, so we get 1/sqrt(n-2) */
    const double tolcorr = std_check / sqrtf((double)count - 2.);

    if ((fabs(mean - 0.5) > tolmean) || (fabs(var - 1. / 12.) > tolvar) ||
        (fabs(correlation) > tolcorr) || (fabs(correlationID) > tolcorr) ||
        (fabs(meanID - 0.5) > tolmean) || (fabs(varID - 1. / 12.) > tolvar) ||
        (fabs(corr_star_sf) > tolcorr) || (fabs(corr_star_se) > tolcorr) ||
        (fabs(corr_star_bh) > tolcorr) || (fabs(corr_sf_se) > tolcorr) ||
        (fabs(corr_sf_bh) > tolcorr) || (fabs(corr_se_bh) > tolcorr) ||
        (fabs(mean_sf - 0.5) > tolmean) || (fabs(mean_se - 0.5) > tolmean) ||
        (fabs(mean_bh - 0.5) > tolmean) || (fabs(var_sf - 1. / 12.) > tolvar) ||
        (fabs(var_se - 1. / 12.) > tolvar) ||
        (fabs(var_bh - 1. / 12.) > tolvar)) {
      message("Test failed!");
      message("Global result:");
      message("Result:    count=%d mean=%f var=%f, correlation=%f", count, mean,
              var, correlation);
      message("Expected:  count=%d mean=%f var=%f, correlation=%f", count, 0.5f,
              1. / 12., 0.);
      message("Max difference:           mean=%f var=%f, correlation=%f",
              tolmean, tolvar, tolcorr);
      message("Difference:               mean=%f var=%f, correlation=%f",
              fabs(mean - 0.5f), fabs(var - 1. / 12.), fabs(correlation));
      message("ID part");
      message(
          "Result:    count=%d mean=%f var=%f"
          " correlation=%f",
          count, meanID, varID, correlationID);
      message(
          "Expected:  count=%d mean=%f var=%f"
          " correlation=%f",
          count, .5f, 1. / 12., 0.);
      message("Max difference:           mean=%f var=%f, correlation=%f",
              tolmean, tolvar, tolcorr);
      message("Difference:               mean=%f var=%f, correlation=%f",
              fabs(meanID - 0.5f), fabs(varID - 1. / 12.), fabs(correlation));
      message("Different physical processes:");
      message(
          "Means:    stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          mean, mean_sf, mean_se, mean_bh);
      message(
          "Expected: stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          .5f, .5f, .5f, .5f);
      message(
          "Max diff: stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          tolmean, tolmean, tolmean, tolmean);
      message(
          "Diff:     stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          fabs(mean - .5f), fabs(mean_sf - .5f), fabs(mean_se - .5f),
          fabs(mean_bh - .5f));
      message(" ");
      message(
          "Var:      stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          var, var_sf, var_se, var_bh);
      message(
          "Expected: stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          1. / 12., 1. / 12., 1 / 12., 1. / 12.);
      message(
          "Max diff: stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          tolvar, tolvar, tolvar, tolvar);
      message(
          "Diff:     stars=%f stellar feedback=%f stellar "
          " enrichment=%f black holes=%f",
          fabs(var - 1. / 12.), fabs(var_sf - 1. / 12.),
          fabs(var_se - 1. / 12.), fabs(var_bh - 1. / 12.));
      message(" ");
      message(
          "Correlation: stars-sf=%f stars-se=%f stars-bh=%f "
          "sf-se=%f sf-bh=%f se-bh=%f",
          corr_star_sf, corr_star_se, corr_star_bh, corr_sf_se, corr_sf_bh,
          corr_se_bh);
      message(
          "Expected:    stars-sf=%f stars-se=%f stars-bh=%f "
          "sf-se=%f sf-bh=%f se-bh=%f",
          0., 0., 0., 0., 0., 0.);
      message(
          "Max diff:    stars-sf=%f stars-se=%f stars-bh=%f "
          "sf-se=%f sf-bh=%f se-bh=%f",
          tolcorr, tolcorr, tolcorr, tolcorr, tolcorr, tolcorr);
      message(
          "Diff:        stars-sf=%f stars-se=%f stars-bh=%f "
          "sf-se=%f sf-bh=%f se-bh=%f",
          fabs(corr_star_sf), fabs(corr_star_se), fabs(corr_star_bh),
          fabs(corr_sf_se), fabs(corr_sf_bh), fabs(corr_se_bh));
      return 1;
    }
  }

  return 0;
}
