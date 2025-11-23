// Auto-generated linear harmonic climate model
// Station name : STRASBOURG-ENTZHEIM
// Station code : 67124001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 7.6403331756591797;
static constexpr double latitude_deg = 48.54949951171875;
static constexpr double delta_utc_solar_h = 0.50935554504394531;
namespace detail {
static constexpr double two_pi = 6.2831853071795864769;
static constexpr double solar_year_days = 365.242189;
static constexpr double omega_annual = two_pi / solar_year_days;
static constexpr double omega_diurnal = two_pi / 24.0;
inline double eval_annual(const double* coeffs, int n_annual, double day){
    double value = coeffs[0];
    for(int k = 1; k <= n_annual; ++k){
        double angle = k * omega_annual * day;
        value += coeffs[2*k - 1] * std::cos(angle);
        value += coeffs[2*k] * std::sin(angle);
    }
    return value;
}
inline double wrap_day(double d){
    while (d >= solar_year_days) d -= solar_year_days;
    while (d < 0.0)   d += solar_year_days;
    return d;
}
inline double wrap_hour(double h){
    while (h >= 24.0) h -= 24.0;
    while (h < 0.0)   h += 24.0;
    return h;
}
} // namespace detail
namespace temperature_model {
static constexpr int n_diurnal = 3;
static constexpr int c0_coeffs_n_annual = 3;
static constexpr double c0_coeffs[] = {
    11.630518845102834, -9.0169417163122016, -2.2770626692190143, 0.022911727235426137, 0.34375688394127463, 0.049366392351824571, -0.086574427153008907
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.4664513117376274, 1.9961347044083024, 0.019136758174910871, 0.59850365845060327, -0.19905601200669373, -0.055503811078381919, 0.053892809269230298
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.0344693451776246, 0.031602630522189437, -0.40667762687014358, 0.20063583801190146, -0.20592067408950843, 0.20171007367172644, 0.11549632216261922
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.51659251346769441, 0.29898275788224749, -0.14929232500035602, -0.31884641453972606, 0.088429681935360646, -0.029554646413761914, 0.03901274637481858
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.013781168429857603, 0.26023602288798064, 0.062156214741882231, 0.11357066877897583, 0.079588480792348371, -0.12776476685951388, -0.0076746976302216703
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.11983024899135239, -0.32004501271537633, 0.043086020747432535, -0.054496629633631971, 0.065426015637783649, 0.046377347929377588, -0.062690222935443407
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.11276502330301716, 0.14666263259093126, 0.004011666032569478, -0.041012575634364112, 0.0023969860282530672, -0.046879502637222667, -0.0052573734727858123
};
inline double evaluate(double day_solar, double hour_solar) {
    double result = detail::eval_annual(c0_coeffs, c0_coeffs_n_annual, day_solar);
    const double ome_d = detail::omega_diurnal;
    result += detail::eval_annual(a1_coeffs, a1_coeffs_n_annual, day_solar) * std::cos(1 * ome_d * hour_solar);
    result += detail::eval_annual(b1_coeffs, b1_coeffs_n_annual, day_solar) * std::sin(1 * ome_d * hour_solar);
    result += detail::eval_annual(a2_coeffs, a2_coeffs_n_annual, day_solar) * std::cos(2 * ome_d * hour_solar);
    result += detail::eval_annual(b2_coeffs, b2_coeffs_n_annual, day_solar) * std::sin(2 * ome_d * hour_solar);
    result += detail::eval_annual(a3_coeffs, a3_coeffs_n_annual, day_solar) * std::cos(3 * ome_d * hour_solar);
    result += detail::eval_annual(b3_coeffs, b3_coeffs_n_annual, day_solar) * std::sin(3 * ome_d * hour_solar);
    return result;
}
} // namespace temperature_model
namespace specific_humidity_model {
static constexpr int n_diurnal = 3;
static constexpr int c0_coeffs_n_annual = 3;
static constexpr double c0_coeffs[] = {
    0.0067930587127140191, -0.0029167275621505605, -0.0014144411595312554, 0.00033321689889889971, 0.00034182715415612114, -2.0271668780245846e-05, 4.350502528914194e-05
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -0.00018020055459292396, -6.3435233346698468e-05, 9.8436698723118824e-05, 9.8455170949358935e-05, 4.1528804106870788e-06, 2.538296590654586e-05, -4.2298345575719862e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -4.7399105712848534e-05, -7.6222117325402855e-05, 1.8822301130500703e-05, 5.3968999044243792e-06, -2.6334196843940134e-05, 1.788859013240576e-05, -7.4068262752256204e-06
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -6.4345949834388247e-05, 0.00016548423913306432, 2.0106442915788485e-05, -4.9643221783930621e-06, -2.9730629134537312e-05, -3.7094714831865813e-05, 3.9836727391837519e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -9.8346142812546155e-05, 0.00010918605855541598, 1.0569933915716571e-05, 2.5665738216245564e-05, -2.095564507921137e-05, -1.069785554308111e-05, 1.0813046248171926e-05
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    1.3072140902738454e-05, 1.8350291667418151e-07, -1.9581077615205618e-05, -5.5721146390391017e-05, 8.42426482565504e-06, 1.466911329948567e-05, 1.1725124015399728e-05
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -2.5718928818810051e-06, 3.4047413895954403e-05, -7.7592270611295876e-06, -1.1049447669802269e-05, 1.1217099455185991e-05, -9.0210309215946717e-06, -1.9347397498457882e-06
};
inline double evaluate(double day_solar, double hour_solar) {
    double result = detail::eval_annual(c0_coeffs, c0_coeffs_n_annual, day_solar);
    const double ome_d = detail::omega_diurnal;
    result += detail::eval_annual(a1_coeffs, a1_coeffs_n_annual, day_solar) * std::cos(1 * ome_d * hour_solar);
    result += detail::eval_annual(b1_coeffs, b1_coeffs_n_annual, day_solar) * std::sin(1 * ome_d * hour_solar);
    result += detail::eval_annual(a2_coeffs, a2_coeffs_n_annual, day_solar) * std::cos(2 * ome_d * hour_solar);
    result += detail::eval_annual(b2_coeffs, b2_coeffs_n_annual, day_solar) * std::sin(2 * ome_d * hour_solar);
    result += detail::eval_annual(a3_coeffs, a3_coeffs_n_annual, day_solar) * std::cos(3 * ome_d * hour_solar);
    result += detail::eval_annual(b3_coeffs, b3_coeffs_n_annual, day_solar) * std::sin(3 * ome_d * hour_solar);
    return result;
}
} // namespace specific_humidity_model
namespace pressure_model {
static constexpr int n_diurnal = 3;
static constexpr int c0_coeffs_n_annual = 3;
static constexpr double c0_coeffs[] = {
    999.0824237561809, 1.1219263742787144, -0.39551158774047257, 0.64603506009202838, 0.56940122112388869, 0.20958431996813401, 0.68719531561001346
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.15832786497652804, -0.13349456845385091, 0.028691663926697114, -0.0090648098275449389, -0.005203894329429846, 0.010757784712490281, 0.010156414904942455
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.32969753134281921, -0.20354469635728559, 0.075585415828834962, -0.04183071509621112, 0.024264237050892392, -0.017635667006553536, -0.0020082991083290977
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.13283093352586839, 0.089244921201677013, 0.041024307665341805, 0.062977214464094239, 0.042067006889176736, -0.048043061208678704, -0.02000776412655067
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.38210143986931477, -0.004552874712299787, -0.02222058453898194, 0.043164114300497373, -0.015876715958359423, 0.017688044330650774, 0.012177488662760359
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.01384515619063365, 0.094998677876742932, -0.035429747036801691, -0.017095135747577924, -0.0094995460784062895, -0.0024613321953398069, 0.0040982800323460748
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.038642475529251319, 0.10088254163286434, 0.0068957672519055907, 0.014794301188151521, 0.012258758457253922, -0.0039008361280434959, 0.0068874858662842479
};
inline double evaluate(double day_solar, double hour_solar) {
    double result = detail::eval_annual(c0_coeffs, c0_coeffs_n_annual, day_solar);
    const double ome_d = detail::omega_diurnal;
    result += detail::eval_annual(a1_coeffs, a1_coeffs_n_annual, day_solar) * std::cos(1 * ome_d * hour_solar);
    result += detail::eval_annual(b1_coeffs, b1_coeffs_n_annual, day_solar) * std::sin(1 * ome_d * hour_solar);
    result += detail::eval_annual(a2_coeffs, a2_coeffs_n_annual, day_solar) * std::cos(2 * ome_d * hour_solar);
    result += detail::eval_annual(b2_coeffs, b2_coeffs_n_annual, day_solar) * std::sin(2 * ome_d * hour_solar);
    result += detail::eval_annual(a3_coeffs, a3_coeffs_n_annual, day_solar) * std::cos(3 * ome_d * hour_solar);
    result += detail::eval_annual(b3_coeffs, b3_coeffs_n_annual, day_solar) * std::sin(3 * ome_d * hour_solar);
    return result;
}
} // namespace pressure_model
namespace stochastic {
// VAR(1) model for residuals: [P, Q, T]
// x_t = x_{t-1} * A + noise
// noise = u * L^T, where u ~ N(0, I)
static constexpr double A[9] = {
    1.1489789580501388, -0.068994595970491013, 0.0026294526510779446,
    -0.0032076078897655541, 0.96310099709522967, 0.062635335388107241,
    -0.0032818887622014792, 0.011425986550274422, 1.0182804812417992,
    -0.15712832225260989, 0.045332551342130545, -0.0085525654027398534,
    0.004775326236789567, -0.013548437163546811, -0.028464162409835739,
    0.0071667764984871592, -0.0019743927610930412, -0.11126692358134821,
    0.0033340380965321409, 0.017888762019912074, 0.015628102232089949,
    -0.001582411357042476, 0.017307943208631771, -0.003598365757941438,
    -0.0029451727552617258, 0.0082104969237816045, 0.047552312740573824
};
static constexpr double L[9] = {
    0.046994270880263332, 0, 0,
    -0.0033709100168280051, 0.20887636472865376, 0,
    -0.042682084446143798, 0.026119065092640895, 0.1932902200437826
};

struct State {
    double p_res = 0.0;
    double q_res = 0.0;
    double t_res = 0.0;
};

// Updates state using VAR(1) process.
// u_p, u_q, u_t: Independent standard normal random numbers provided by the caller.
inline void update_state(State& s, double u_p, double u_q, double u_t) {
    // Map state and input to arrays
    double x[3] = {s.p_res, s.q_res, s.t_res};
    double u[3] = {u_p, u_q, u_t};

    // Compute noise n = u * L^T
    // n[j] = sum_k u[k] * L[j][k] (since L^T at [k][j] is L[j][k])
    double n[3] = {0.0, 0.0, 0.0};
    for(int j=0; j<3; ++j) {
        for(int k=0; k<3; ++k) {
             n[j] += u[k] * L[j*3 + k];
        }
    }

    // Compute transition v = x * A
    // v[j] = sum_k x[k] * A[k][j] (assuming row vector x)
    double v[3] = {0.0, 0.0, 0.0};
    for(int j=0; j<3; ++j) {
        for(int k=0; k<3; ++k) {
             v[j] += x[k] * A[k*3 + j];
        }
    }

    s.p_res = v[0] + n[0];
    s.q_res = v[1] + n[1];
    s.t_res = v[2] + n[2];
}

} // namespace stochastic
inline double predict_temperature(double day_utc, double hour_utc){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    return temperature_model::evaluate(day_solar, hour_solar);
}
inline double predict_specific_humidity(double day_utc, double hour_utc){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    return specific_humidity_model::evaluate(day_solar, hour_solar);
}
inline double predict_pressure(double day_utc, double hour_utc){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    return pressure_model::evaluate(day_solar, hour_solar);
}
inline void predict(double day_utc, double hour_utc, double& temperature_c, double& specific_humidity_kg_kg, double& pressure_hpa){
    double hour_solar = detail::wrap_hour(hour_utc + delta_utc_solar_h);
    double day_solar  = detail::wrap_day(day_utc + (delta_utc_solar_h / 24.0));
    temperature_c = temperature_model::evaluate(day_solar, hour_solar);
    specific_humidity_kg_kg = specific_humidity_model::evaluate(day_solar, hour_solar);
    pressure_hpa = pressure_model::evaluate(day_solar, hour_solar);
}
inline void predict_stochastic(double day_utc, double hour_utc, const stochastic::State& state, double& temperature_c, double& specific_humidity_kg_kg, double& pressure_hpa){
    predict(day_utc, hour_utc, temperature_c, specific_humidity_kg_kg, pressure_hpa);
    temperature_c += state.t_res;
    specific_humidity_kg_kg += state.q_res;
    pressure_hpa += state.p_res;
}
} // namespace harmoclimat