
# package imports
import numpy as np

# 4% enriched at 50 MWd/kgU
decay_data = { "U-235":
                       { "power_fraction": 0.235,
                         "energy_per_fission" : 202.2,
                         "alpha" : np.array([ 5.2800E-04,
                                     6.8588E-01,
                                     4.0752E-01,
                                     2.1937E-01,
                                     5.7701E-02,
                                     2.2530E-02,
                                     3.3392E-03,
                                     9.3667E-04,
                                     8.0899E-04,
                                     1.9572E-04,
                                     3.2609E-05,
                                     7.5827E-06,
                                     2.5189E-06,
                                     4.9836E-07,
                                     1.8523E-07,
                                     2.6592E-08,
                                     2.2356E-09,
                                     8.9582E-12,
                                     8.5968E-11,
                                     2.1072E-14,
                                     7.1219E-16,
                                     8.1126E-17,
                                     9.4678E-17]),
                         "lambda" : np.array([ 2.7216E+00,
                                      1.0256E+00,
                                      3.1419E-01,
                                      1.1788E-01,
                                      3.4365E-02,
                                      1.1762E-02,
                                      3.6065E-03,
                                      1.3963E-03,
                                      6.2608E-04,
                                      1.8924E-04,
                                      5.5074E-05,
                                      2.0971E-05,
                                      9.9940E-06,
                                      2.5401E-06,
                                      6.6332E-07,
                                      1.2281E-07,
                                      2.7163E-08,
                                      3.2955E-09,
                                      7.4225E-10,
                                      2.4681E-10,
                                      1.5596E-13,
                                      2.2573E-14,
                                      2.0503E-14])},
               "Pu-239":
                       { "power_fraction" : 0.512,
                         "energy_per_fission" : 205.5,
                         "alpha" : np.array([ 3.9368E-01,
                                     7.4588E-01,
                                     1.2169E+00,
                                     5.2820E-01,
                                     1.4805E-01,
                                     4.5980E-02,
                                     1.0406E-02,
                                     1.6991E-03,
                                     6.9102E-04,
                                     1.4736E-04,
                                     2.4049E-05,
                                     6.9288E-06,
                                     6.4927E-07,
                                     4.3556E-07,
                                     1.6020E-07,
                                     2.3089E-08,
                                     2.5481E-09,
                                     3.5071E-11,
                                     6.3399E-11,
                                     4.1599E-14,
                                     5.3295E-16,
                                     1.6695E-18,
                                     4.1058E-16]),
                         "lambda" : np.array([ 4.3427E+00,
                                      1.7114E+00,
                                      6.0572E-01,
                                      1.9429E-01,
                                      6.9788E-02,
                                      1.8809E-02,
                                      6.1265E-03,
                                      1.3799E-03,
                                      5.2799E-04,
                                      1.6145E-04,
                                      4.8419E-05,
                                      1.5644E-05,
                                      5.3610E-06,
                                      2.1689E-06,
                                      6.3343E-07,
                                      1.2879E-07,
                                      2.5604E-08,
                                      9.1544E-09,
                                      7.3940E-10,
                                      2.4731E-10,
                                      1.9594E-13,
                                      6.4303E-14,
                                      6.4229E-14])},
               "U-238":
                       { "power_fraction" : 0.091,
                         "energy_per_fission" : 205.5,
                         "alpha" : np.array([ 3.9368E-01,
                                     7.4588E-01,
                                     1.2169E+00,
                                     5.2820E-01,
                                     1.4805E-01,
                                     4.5980E-02,
                                     1.0406E-02,
                                     1.6991E-03,
                                     6.9102E-04,
                                     1.4736E-04,
                                     2.4049E-05,
                                     6.9288E-06,
                                     6.4927E-07,
                                     4.3556E-07,
                                     1.6020E-07,
                                     2.3089E-08,
                                     2.5481E-09,
                                     3.5071E-11,
                                     6.3399E-11,
                                     4.1599E-14,
                                     5.3295E-16,
                                     1.6695E-18,
                                     4.1058E-16]),
                         "lambda" : np.array([ 4.3427E+00,
                                      1.7114E+00,
                                      6.0572E-01,
                                      1.9429E-01,
                                      6.9788E-02,
                                      1.8809E-02,
                                      6.1265E-03,
                                      1.3799E-03,
                                      5.2799E-04,
                                      1.6145E-04,
                                      4.8419E-05,
                                      1.5644E-05,
                                      5.3610E-06,
                                      2.1689E-06,
                                      6.3343E-07,
                                      1.2879E-07,
                                      2.5604E-08,
                                      9.1544E-09,
                                      7.3940E-10,
                                      2.4731E-10,
                                      1.9594E-13,
                                      6.4303E-14,
                                      6.4229E-14])},
               "Pu-241":
                       { "power_fraction" : 0.162,
                         "energy_per_fission" : 205.5,
                         "alpha" : np.array([ 3.9368E-01,
                                     7.4588E-01,
                                     1.2169E+00,
                                     5.2820E-01,
                                     1.4805E-01,
                                     4.5980E-02,
                                     1.0406E-02,
                                     1.6991E-03,
                                     6.9102E-04,
                                     1.4736E-04,
                                     2.4049E-05,
                                     6.9288E-06,
                                     6.4927E-07,
                                     4.3556E-07,
                                     1.6020E-07,
                                     2.3089E-08,
                                     2.5481E-09,
                                     3.5071E-11,
                                     6.3399E-11,
                                     4.1599E-14,
                                     5.3295E-16,
                                     1.6695E-18,
                                     4.1058E-16]),
                         "lambda" : np.array([ 4.3427E+00,
                                      1.7114E+00,
                                      6.0572E-01,
                                      1.9429E-01,
                                      6.9788E-02,
                                      1.8809E-02,
                                      6.1265E-03,
                                      1.3799E-03,
                                      5.2799E-04,
                                      1.6145E-04,
                                      4.8419E-05,
                                      1.5644E-05,
                                      5.3610E-06,
                                      2.1689E-06,
                                      6.3343E-07,
                                      1.2879E-07,
                                      2.5604E-08,
                                      9.1544E-09,
                                      7.3940E-10,
                                      2.4731E-10,
                                      1.9594E-13,
                                      6.4303E-14,
                                      6.4229E-14])}}
