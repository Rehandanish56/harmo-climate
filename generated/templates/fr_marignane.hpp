// Auto-generated linear harmonic climate model
// Station name : MARIGNANE
// Station code : 13054001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 5.2160000801086426;
static constexpr double latitude_deg = 43.437667846679688;
static constexpr double delta_utc_solar_h = 0.34773333867390954;
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
    15.976921618622223, -8.5613983703078951, -3.197460083525514, 0.064455888487419094, 0.6405240341717795, -0.35967327582975217, -0.13229505400198818
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.5558472342713205, 1.1301622180809039, 0.011639389758802712, 0.1604754824620635, -0.16085338785775322, -0.043556446367857458, 0.0099817718842041426
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -0.83716920120529725, -0.5443527391791203, -0.24331951450517672, -0.033882816435339705, -0.31935462586044522, 0.1695069413432255, 0.054512393672620542
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.65609671480883613, 0.38446608430134943, -0.075746168228386854, -0.046948699076037735, 0.10602569344385217, -0.04266525475383217, -0.041138258818060193
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.2608872596039753, 0.47915191246786837, 0.10312840194564848, 0.15930817736252889, 0.050130979706688089, -0.073641778598143678, 0.014909540291882016
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.06210684266096702, -0.22156872531502733, 0.023269203219176176, -0.097101770082856004, 0.012006898145627691, 0.00060430444232454391, -0.0082835018315489239
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.080881693915711705, 0.18102168938086283, -0.03564027179882872, -0.058724842783423319, 0.013847774824410339, -0.050111838404602281, -0.016866422501775471
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
    0.0075912077323257218, -0.0026947116466066751, -0.001621417447463978, 5.7681695759935985e-05, 0.00021864774757915248, -0.00014086746762538441, 6.706258431641634e-05
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -0.00014107820225910627, -0.00019392114481041312, 2.3275445025389654e-06, 0.00010785515326082917, 7.8728724549001982e-05, -5.0145791294626327e-06, -2.9227155471867144e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -0.00013441672866548876, 5.0773201203074643e-06, 5.2850273628071587e-05, -2.3553372369081156e-05, -7.1525182924090768e-05, 5.4961663956482404e-06, 1.4804179783523134e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -7.5702837962837629e-05, 0.00017243808103166019, 2.9666906830506811e-05, 4.9161305075963638e-06, -3.6739066412980104e-05, -1.8315317107108201e-05, 3.8328027137733012e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -8.726548854274536e-05, 6.2882673730573797e-05, -4.1580884054441979e-06, 1.7173250424570115e-05, -8.8146682522911939e-06, -2.8012527635433158e-06, 7.5746165866309795e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    2.5004208934600869e-05, 1.0126299806409624e-05, -2.51018804145144e-05, -4.4355958651116732e-05, 1.4327113516133925e-05, 2.152684134050188e-06, -2.9687961798338829e-07
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -4.1155040324719033e-06, 5.6339547825705959e-05, -2.7219303466289045e-06, -5.7968409642453559e-06, 8.2117679676296862e-06, -1.1626830191867437e-05, -5.5862124593631963e-06
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
    1012.1520112242902, 1.6639678217220875, -0.27368016380893351, 0.75589559238288362, 0.38262622917102235, 0.2497420883605731, 0.69142808211115048
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.22969398395830695, -0.053700172637602295, 0.012127698458284, -0.012816517326044488, 0.0027516643533722989, -0.0024155820599033361, -0.0076742412337015191
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.19412381034557916, -0.066949820259036619, 0.021661167020933268, 0.020788371064179922, 0.078692987366580938, -0.022216844657134509, -0.015630889220692914
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.14580663007053635, 0.12000152286964669, 0.06158887088145483, 0.066488447458881411, 0.045659142336767253, -0.047527965065948058, -0.012702484988989401
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.51049450376878758, -0.066851787383041428, -0.03428834605836701, 0.037258031621822077, -0.017553190045870068, 0.014074441057477887, 0.011790775110123075
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.011749852634012322, 0.12437781474967947, -0.037079286514683923, -0.018930085831049416, -0.0090095481925165921, -0.0061035617544220981, 0.0084578228795940662
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.055021469385746008, 0.11941814163025, 0.01184441858715057, 0.026625874195932885, 0.016660383419900348, -0.0030535730185875117, 0.006651310860486149
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
    1.0772108995021719, -0.05287646605152091, 0.002021645452802616,
    -0.0066721728800557059, 0.93433801692220242, 0.038779217135915402,
    -0.0011598328080020454, 0.019087479929335877, 0.98246910602242588,
    -0.045046698576412005, 0.053441617067625735, -0.0071802141376450176,
    0.0039993706795369193, 0.015728361500801712, -0.016658032685972872,
    0.0052853509455702673, 0.005865337564117265, -0.08426841784947986,
    -0.038334349009529617, 0.0020862688749396785, 0.018934394851393411,
    0.00078588425054546994, 0.0067971715875004834, 0.010229920012323207,
    -0.0022831819411089542, -0.0033628619989132073, 0.047596015807329736
};
static constexpr double L[9] = {
    0.066922143109720944, 0, 0,
    0.010321177778319056, 0.25889872031519273, 0,
    -0.043818187074847685, 0.0094936784361589841, 0.25636497720628332
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