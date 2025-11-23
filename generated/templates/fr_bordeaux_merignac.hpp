// Auto-generated linear harmonic climate model
// Station name : BORDEAUX-MERIGNAC
// Station code : 33281001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = -0.69133299589157104;
static constexpr double latitude_deg = 44.830665588378906;
static constexpr double delta_utc_solar_h = -0.046088866392771419;
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
    14.091954666830921, -7.2344516508833729, -2.5975360719274621, -0.21383915029157227, 0.61018563101483736, -0.11381710749095345, 0.017048700134991483
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.5664500746977303, 1.5322911406076078, 0.16890533846691519, 0.43186745300664969, -0.229667740573992, -0.042500649358237746, 0.010539215953296759
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -0.91856981456273257, -0.31301483110304651, -0.27352309889779997, 0.053179589773275637, -0.2443072499675509, 0.18863123759905032, 0.067935978724252777
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.62210158726244746, 0.36170036680480844, -0.12915009843645908, -0.16736436837747373, 0.044007308183882793, -0.054638232403163178, -0.011402241623500481
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.067458198144992293, 0.28439883627181872, 0.10250824062948315, 0.18830620519471461, 0.053017666239912338, -0.075819891399385525, -0.00069441491375159085
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.13865602238372851, -0.30830680932938631, 0.016723539273476966, -0.10965329056190098, 0.052101773844991602, 0.033196405631412504, -0.01007079523279718
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.10984556414116636, 0.18648894174615832, 0.0074967980117907219, -0.052295324924943203, 0.0090971643417298402, -0.079684628238395644, -0.007706169510713309
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
    0.0076808282907577501, -0.0023896052396574709, -0.0014409766913201212, 0.00020008100914717229, 0.00023418660530284628, -0.00010958204125407822, 0.00010887146376858432
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -2.7952359795197422e-05, -0.00023642081831479, 2.8871735647062328e-05, 2.6595610012323663e-05, 9.1545925203138269e-05, 5.7711191051293882e-05, -5.2710880219493432e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.8773892935316445e-05, -0.00014224417851758386, 7.2494442962968954e-07, -1.8961460642700359e-05, 3.1051949932983143e-05, 1.8509071241619363e-05, -3.6218511722187733e-06
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -7.9095589402439833e-05, 0.00014746713341592119, 2.1096915623021741e-05, 4.0661925717598106e-05, -2.2589787110005856e-05, -3.4048887073531448e-05, -8.6890049469545156e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -8.7709790040784512e-05, 4.411331425730938e-05, 8.343874134491929e-06, 4.5139291546019073e-05, -2.278586662869658e-05, 9.195817140995816e-06, -5.210509345335214e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    1.4367168033024098e-05, 2.4543338105129348e-05, -2.5043640591647111e-05, -5.4446766983561142e-05, -3.1999459788070989e-08, -3.1817426054056316e-06, 1.3174709419141795e-05
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -5.498845697001555e-06, 6.7793537415721482e-05, 3.1898197535262822e-06, -3.7955646306291599e-06, -6.8806809543801335e-07, -2.0363947541256718e-05, 1.7486044641318267e-06
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
    1011.7632789588021, 1.1606943901357145, -0.076521241810292254, 1.4101296454181129, 0.91827266802542795, 0.36643331026827325, 0.6551672144841969
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.10606776502440471, -0.10647817656794435, 0.031089860180092282, -0.006558830451799909, 0.029471092539713205, 0.0085197599230108934, -0.010660241860948404
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.10064850604435159, -0.096556369006179255, 0.018579042943014632, -0.026768101727461736, 0.044548230364737117, -0.00021614711152414691, -0.011402581923263277
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.18876193591559823, 0.10943626606786536, 0.056273362975654717, 0.064160067693367712, 0.05610556541687288, -0.053972048283000777, -0.013242912236813582
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.46585922165057519, -0.028904231234860925, -0.028413389719962833, 0.022566442386888, -0.018140765001737037, 0.010758737632605604, 0.015392895688649629
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.013257859894161561, 0.13218145879553739, -0.033464833121529272, -0.021558666576108037, -0.0035919144899792112, -0.0087972267658650068, 0.0051924840530567541
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.049552506309794714, 0.10673548634871598, 0.012186121723048705, 0.022063714775340429, 0.013853057713785267, 0.0040306549897806436, 0.0073475884028625918
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
    1.1473534303930943, -0.092835321000514448, -0.05347389904819573,
    7.828477901213049e-05, 0.97267862972250752, 0.048503167357045032,
    -0.001813969813479293, 0.0026499618622618376, 1.0150836158932739,
    -0.15140826329168949, 0.1054143392130199, 0.065985230574981207,
    0.0020589252880991327, -0.022991694979389854, -0.02500066797990335,
    0.0061678515486029203, 0.019423326453347933, -0.1017057941019925,
    -0.00091244366639601365, -0.020323851145995886, -0.0048757755037848044,
    -0.0018047396980928904, 0.011966053958390704, 0.00037226433511945364,
    -0.003660559865074544, 0.0012472487584677117, 0.045725322060144895
};
static constexpr double L[9] = {
    0.048791629962126555, 0, 0,
    0.0018236986811009496, 0.21460134585118329, 0,
    -0.042818666089719862, 0.011659531514045992, 0.20099786975667042
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