def cusp_forms_dim(p, k):
    weight_par = (-1) ** k
    weight_over_three = [0, 1, -1][k % 3]
    weight_six = [2, 1, -1, -2, -1, 1][k % 6]
    kronecker_4 = kronecker_symbol(-1, p)
    kronecker_3 = kronecker_symbol(-3, p)
    alpha_0 = (p + 1) * (p^2 + 1) * (2 * k - 2) * (2 * k - 3) * (2 * k - 4) / 69120
    alpha_1 = (p + 2 + kronecker_4) * weight_par / 128
    alpha_2_3 = -weight_over_three / 108
    if p > 3:
        alpha_2_3 *= (p + 2 + kronecker_3)
    if p % 8 == 7:
        alpha_4_5_6 = 0
    elif p % 8 == 3 or p % 8 == 5:
        alpha_4_5_6 = 1 / 4
    else:
        alpha_4_5_6 = 1 / 2
    if k % 4 == 3:
        alpha_4_5_6 *= -1
    elif k % 4 == 1 or k % 4 == 2:
        alpha_4_5_6 = 0
    alpha_7_8 = (1 + kronecker_3)^2 * weight_par / 27
    alpha_9_10_11_12 = (1 + kronecker_3)^2 * weight_six / 27
    alpha_13_14 = weight_over_three * [0, 4, 0, 1, 0, 2, 0, 2, 0, 0, 0, 0][p % 12] / 12
    alpha_15_16_17_18 =  [1, 0, 0, -1, 0][k % 5] * [1, 4, 0, 0, 0][p % 5] / 5
    alpha_19_20_21_22 = [1, 0, 0, -1, -1, -1, -1, 0, 0, 1, 1, 1][k % 12] * (1 + kronecker_4) * (1 + kronecker_3) / 12
    beta_1_2 = [2*k - 3, -k + 1, -k + 2][k % 3] * (p + 1) * (1 + kronecker_3) / 216
    beta_3_4 = [-1, -k + 1, -k + 2, 1, k - 1, k - 2][k % 6] * (p + 1) * (1 + kronecker_3) / 72
    beta_5_6 = [k - 2, -k + 1, -k + 2, k - 1][k % 4] * (p + 1) * (1 + kronecker_4) / 96
    gamma_1_2 = 5 * (2 * k - 3) * (p + 2 + kronecker_4) / 384
    gamma_3 = (2 * k - 3) * [7, p + 2 + kronecker_3][p > 3] / 54
    delta_1_2 = 7 * weight_par * (2 * k - 2) * (2 * k - 4) * (p + 1)^2 / 4608
    beta_hat_1_2 = [0, 1, 1, 0, -1, -1][k % 6] * (1 + kronecker_3) / 6
    beta_hat_3_4 = [-2, 1, 1][k % 3] * (1 + kronecker_3) / 18
    if p == 3:
        beta_hat_3_4 += [-1, 1, 0][k % 3] / 9
    beta_hat_5_6 = [-1, 1, 0][k % 3] * (1 + kronecker_3) * [1, 2][p > 3] / 9
    beta_hat_7_8 = [-1, 1, 1, -1][k % 4] * (1 + kronecker_4) / 8
    beta_hat_9_10 =  [-1, 1, 1, -1][k % 4] * (1 + kronecker_4) / 8
    delta_hat_hat_1_2 = weight_par / 4
    delta_hat_hat_3_4 = weight_par * [0, 1, 0, 2][p % 4] / 8
    delta_hat_1_2 = weight_par * (3 - 2 * k) * (p + 1) / 24
    eps_1_2_3 = 1/12
    eps_4 = (3 - 2 * k) * (p + 1) / 144
    gamma_hat_1_2_3_4 = (3 + kronecker_4) / (-8)
    gamma_hat_5_6_7 = (3 + kronecker_3) / (-6)
    return alpha_0 + alpha_1 + alpha_2_3 + alpha_4_5_6 + alpha_7_8 + alpha_9_10_11_12 + alpha_13_14 + alpha_15_16_17_18 + alpha_19_20_21_22 + beta_1_2 + beta_3_4 + beta_5_6 + gamma_1_2 + gamma_3 + delta_1_2 + beta_hat_1_2 + beta_hat_3_4 + beta_hat_5_6 + beta_hat_7_8 + beta_hat_9_10 + delta_hat_hat_1_2 + delta_hat_hat_3_4 + delta_hat_1_2 + eps_1_2_3 + eps_4 + gamma_hat_1_2_3_4 + gamma_hat_5_6_7

def modular_forms_dim(p, k):
    if k % 2:
        return cusp_forms_dim(p, k)
    s = CuspForms(Gamma0(p), k).dimension()
    return cusp_forms_dim(p, k) + s * 2 + 3