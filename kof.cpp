for (int m = 1; m < M; m++)
        {
            mat._main.push_back (
                                 H_n_1[m] / t
                                 + 2 * Mu / (h * h)
            );

            mat._down.push_back (
                                 -Mu / (h * h)
                                 - H_n_1[m - 1] * V_n[m - 1] / (3 * h)
                                 - H_n_1[m] * V_n[m] / (3 * h)
            );

            mat._up.push_back(
                               -Mu / (h * h)
                               + H_n_1[m + 1] * V_n[m + 1] / (3 * h)
                               + H_n_1[m] * V_n[m] / (3 * h)
            );

            b.push_back(
                        H_n_1[m] * f(i, m)
                        + H_n[m] * V_n[m] / t
                        - (1/3. * V_n[m] * V_n[m] + C) * (H_n_1[m + 1] - H_n_1[m - 1]) / (2 * h)
                        
            );
        }