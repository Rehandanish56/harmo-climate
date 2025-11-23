// Auto-generated linear harmonic climate model
// Station name : NANTES-BOUGUENAIS
// Station code : 44020001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = -1.6088329553604126;
static constexpr double latitude_deg = 47.150001525878906;
static constexpr double delta_utc_solar_h = -0.10725553035736081;
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
    12.728277046551483, -6.4904256123643345, -2.4704752648870572, -0.0071601456142067988, 0.47963271087063591, -0.076672208374456152, 0.047973560721237554
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.2696716373138011, 1.6284541172560312, 0.12920379953930328, 0.38602125292938444, -0.18730689877185194, -0.093625024042979454, 0.024068328873576128
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -0.63629561368549437, -0.25278473003059798, -0.23544111324280262, 0.027873860617909994, -0.17785842645861738, 0.13625216900449091, 0.097092715533671548
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.56778466530364446, 0.32468602865077029, -0.15313021721722844, -0.1824245373936399, 0.079898124848041555, -0.056416633881070613, -0.014173861042839069
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.1469316021155474, 0.24244043081615255, 0.099674468815603573, 0.18444958834074954, 0.052997784613093679, -0.070135449622297577, -0.017733695683667297
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.09757408492927469, -0.29437013468882439, 0.030111880970079885, -0.1047067968873244, 0.048718716237224687, 0.038354248840717343, -0.019382100465055913
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.096999725541299639, 0.20176826555812591, -0.027042710730379969, -0.044643798621603621, -0.0013088192446828073, -0.070524108805052332, -0.0073762021878056167
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
    0.0073525296056313446, -0.0021076406216363535, -0.0013180328977171162, 0.00017622602476154018, 0.0002256610407764581, -6.8101986354861642e-05, 7.368175120610958e-05
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -0.00020924107033484278, -4.9047739422899835e-05, 2.0291524676504362e-06, 3.7381440513833622e-05, 7.8718844337620831e-05, 3.1606595471637103e-05, -3.3714740599358867e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -3.5041020085214916e-06, -0.00010604935484995132, -2.1039846207261736e-05, -1.3924895038575729e-07, 1.8583328170159353e-05, 7.0860438752767337e-06, 3.8426218568489887e-06
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -4.6933409827199466e-05, 0.00015747126367192578, 3.3856431445253943e-05, 1.8077934965991283e-05, -3.7288640081608911e-05, -3.5335621141400223e-05, -1.7628199156743553e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -7.511785590059197e-05, 2.4880270819303149e-05, 1.4208022697681542e-05, 4.673015054653205e-05, -5.2159686019266028e-06, 2.0874018447154934e-06, -2.5888318731284699e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    1.392043329488807e-05, 2.2091534918114715e-06, -2.4179255968963558e-05, -5.2889282287667452e-05, 9.3953767976720574e-06, -9.4203886509610605e-07, 1.0831458735855775e-05
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -1.9029314493679949e-05, 7.4333552690269549e-05, 4.2850731565757413e-06, 9.4142638023928255e-07, -7.8962263416688778e-06, -2.0878192324694313e-05, -5.8608836583681195e-06
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
    1014.0554798882262, 0.5820689155897929, 0.082911870944819679, 1.1171534558285381, 0.96482124083588716, 0.45172095792553413, 0.66855427704772763
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.082287097855057015, -0.062079590574974065, 0.026566755793673364, 0.0059205293306923987, 0.013552195641468499, 0.014986270114466895, -0.0028548831348460062
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.08315943799329685, -0.11598592163017594, 0.006745443545817356, -0.011209724445439955, 0.017535343800319944, 0.013407345113429177, -0.013176633237257413
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.14733644655042444, 0.083693699024243298, 0.055137875170536213, 0.063669028794528751, 0.045880190084852022, -0.045848588405960448, -0.018469719578958693
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.40972920785411238, -0.045366432215387317, -0.025454246976593677, 0.021184301157121586, -0.020699257156985656, 0.0085289078989015093, 0.0058758881412002886
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.010086503990821597, 0.1113414341339843, -0.034795257723319178, -0.020657783375799901, -0.0024717278286558973, -0.0092888663583202344, 0.010424329960640856
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.039012074764222415, 0.10733115825766221, 0.0079448444744933486, 0.020496735343457854, 0.013757508038436462, -0.0014788428457152428, 0.003187717759048339
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
    1.1664552479830077, -0.16121787632273124, -0.092500059234021326,
    -0.00082641419270122751, 0.97201627116778222, 0.048905776969712204,
    0.0010997210802492219, -0.0076571124672158777, 0.99058614807332313,
    -0.18107406512610705, 0.21439573630020115, 0.14453416127898425,
    0.0020902064268671118, -0.016626989045893159, -0.013090515766117382,
    0.0022028078011407954, 0.019863367130644252, -0.086789313396311726,
    0.01059638167906031, -0.059630099315456486, -0.041087215369639926,
    -0.00056556064156093253, 0.012493846958941374, -0.004504242645907785,
    -0.003080806501557393, -0.00047342439122543406, 0.046341672329614926
};
static constexpr double L[9] = {
    0.040568002909863912, 0, 0,
    -0.012569331832119647, 0.2186413001447913, 0,
    -0.039932988987216214, 0.052926589666600583, 0.21980799981519464
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