// Auto-generated linear harmonic climate model
// Station name : BOURGES
// Station code : 18033001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 2.3598330020904541;
static constexpr double latitude_deg = 47.059165954589844;
static constexpr double delta_utc_solar_h = 0.15732219815254211;
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
    12.340359813713814, -7.9840292317407835, -2.6092246842536855, 0.0020406793803234269, 0.64036004587844719, -0.12088322905697929, -0.10676041525179968
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.3280905971402199, 1.7340646154812642, 0.062746997428949652, 0.53135756667232992, -0.15142111415419882, -0.055522899832822166, 0.0042939937519765447
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.1915294914645549, 0.11769209371261294, -0.2166770054350374, 0.11951274867407305, -0.32889548122046364, 0.19055693766624751, 0.14091889893412526
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.53367003254708045, 0.3299129885311613, -0.11299893961904016, -0.27469719351346689, 0.060547857514974428, -0.05571261902881116, 0.014894392028617255
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.12117777191870953, 0.33403147282589768, 0.06485561435647097, 0.15384299039445751, 0.062701132012918048, -0.093718616323040196, -0.012282006580291457
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.13340701132360894, -0.33378648648198339, 0.018111748070319969, -0.083708636402585523, 0.047633945640728814, 0.048772810727360609, -0.02919669441838001
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.097975338234572981, 0.16607449027423934, -0.021520596415070741, -0.06264645106324683, 0.019260177792027978, -0.067698004088251537, -0.0042639948185703809
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
    0.0069939521025643983, -0.0023714750349956358, -0.0012722720895884576, 0.00021816022543698726, 0.00012262807104116368, -7.9940493718415594e-05, 0.00015314897003319185
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -0.00013115710550463009, -0.00014035985700408432, 6.1597431470467143e-06, 9.8308063666869948e-05, 9.2487513086479185e-05, 1.4316164394291843e-05, -5.0902607603140491e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -4.3915343000838671e-05, -0.00011125809103588861, -2.542832853025811e-05, 1.1428484061795104e-05, 6.0667543937161684e-05, 3.041577719518654e-05, -4.3090092680687407e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -6.4353976025069561e-05, 0.00017253806544363832, 1.2814570965262049e-05, -2.0349139816186936e-06, -1.9126093326587933e-05, -3.481798889616745e-05, -1.1980576016551576e-05
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -9.124076051098157e-05, 7.5370365464592839e-05, 1.8336455049007561e-05, 5.1235933222299073e-05, -1.9407839474620293e-05, -9.4598653717837483e-06, -1.6950632258372751e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    2.4833427699714072e-05, -1.1493362305739041e-05, -1.7403131141292821e-05, -5.6068876840298624e-05, 5.021565624807691e-06, 1.0354444563596218e-05, 1.0281163925456041e-05
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -9.722181894505486e-06, 5.7436581251018959e-05, -3.5883578170117407e-06, -1.3342679886325975e-05, 5.0365519228834131e-06, -2.2520253253804813e-05, -3.0510278346354022e-06
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
    997.61421650600096, 0.81433658525371366, -0.22121132347864678, 1.0154651620079509, 0.78945953196389007, 0.36714269953549233, 0.67900295280635947
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.11517676100293973, -0.068148539039582332, 0.023196602005184589, 0.0020935350775873108, 0.0055714757886034475, -0.0042669131661350473, 0.0010725267990831459
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.23058386102714462, -0.20779485607484996, 0.023748515775611601, -0.015773369830955133, 0.052646508748258727, -0.015103712027374355, -0.036832356427670729
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.13589111576220833, 0.093920010682813457, 0.041826251477158301, 0.064074774259933451, 0.046015778352966524, -0.041033995542662433, -0.019810963366222949
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.41531781980729832, -0.022278898581923643, -0.016338607605656346, 0.042401602548647464, -0.015696013638781519, 0.0095609822134987265, 0.017600529837321376
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.013853218971503128, 0.10813708271191537, -0.035091801171902477, -0.017116559458609824, -0.0092405457121331816, -0.0031323007271190834, 0.0076453235308270411
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.044856833887575616, 0.10283571613581605, 0.0037964056460657934, 0.019369258241526233, 0.010493870271583957, -0.0013456325058235699, 0.0044893605508164059
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
    1.151669003573949, -0.095364537958633158, -0.021124437589694005,
    -0.0016385917589136906, 0.98808400353900261, 0.052497480729355533,
    -0.0015004744659219424, -0.0027832165370486259, 1.0141031530307394,
    -0.15797632427282102, 0.10446287218462713, 0.022606225660004961,
    0.0039612968877545396, -0.030099802507772472, -0.024079700859695775,
    0.0042794757501498845, 0.018155354842754781, -0.097874137614764198,
    0.0016206230840838318, -0.015405568895721533, 0.0086413258644063504,
    -0.0019070856665606775, 0.010564793599866507, -0.0065331503089491073,
    -0.0023099865974608941, 0.00084374803109932284, 0.046722745985654637
};
static constexpr double L[9] = {
    0.046104230069742784, 0, 0,
    -0.0030364641025691311, 0.19486502229674393, 0,
    -0.040164654583411973, 0.027992212632303011, 0.18981843921244271
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