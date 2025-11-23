// Auto-generated linear harmonic climate model
// Station name : LYON-BRON
// Station code : 69029001
#pragma once
#include <cmath>
namespace harmoclimat {
static constexpr double longitude_deg = 4.9491667747497559;
static constexpr double latitude_deg = 45.721332550048828;
static constexpr double delta_utc_solar_h = 0.32994445164998376;
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
    13.283013490565564, -9.086519622100008, -2.7143910089665275, -0.10370375253688795, 0.65057960052047192, -0.24632173655900413, -0.18523321713124122
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -3.2631302923332011, 1.6141749664552851, 0.036834459604233268, 0.40275484435004816, -0.22719807155949631, -0.027614933365050925, 0.026026839506185308
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -1.1982258228628195, 0.034025518754278603, -0.2689316135655051, 0.14877658307534652, -0.29375772129197386, 0.19793935549210273, 0.1279829873908179
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    0.56932375068429475, 0.27626286737261507, -0.090763892498842272, -0.20483364218715075, 0.067363322500066436, -0.048666062797097627, 0.010350274838791295
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.085622644100533829, 0.33447779613043788, 0.060461976186481414, 0.12234404601892508, 0.062040634773905541, -0.069523648967103779, -0.0093520286357408038
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.106830673584002, -0.26113686341397419, 0.040568114389684835, -0.08404857627411709, 0.055603274184998483, 0.033232331877305975, -0.029041705246676271
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    -0.11834535693534642, 0.16922002957547727, -0.015252765413806825, -0.046199838350084711, 0.020879908687818697, -0.055244069751943362, -0.010059237955686308
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
    0.0068835816433089137, -0.0026366060977896875, -0.0014094113492455951, 0.0001543298998427312, 0.00018361564873233368, -8.0883802419683736e-05, 0.00013171745739984638
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    -1.5401769469129419e-05, -0.00019532296311246388, -1.3267152157274986e-05, 5.3999724795428972e-05, 9.5561122759305476e-05, 1.6337171074995095e-05, -5.2528693177806851e-05
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    -4.0214703058553525e-05, -8.3617204365972478e-05, -1.3027426430393606e-05, 6.8055690563286557e-06, 4.4896770796295516e-05, 1.9861371505555945e-05, -3.2873190055556059e-05
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -5.3091200375911112e-05, 0.00011664600220418934, 2.2649491351288943e-05, 1.4026359019832038e-05, -1.6861085226651931e-05, -1.8500858226132548e-05, -7.3892487559601282e-06
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -8.321344619094776e-05, 6.2661490008277974e-05, 7.4717204737915032e-06, 3.4897337183320911e-05, -1.8536054293554221e-05, -5.7751335574290303e-06, 2.1688898669045863e-06
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    1.563093136285931e-05, -3.0495174797951579e-06, -1.2882582723153151e-05, -2.2403910674672017e-05, 4.4791237476051338e-06, -1.5911704805683714e-06, 8.5423738530742141e-06
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    3.5110664493651384e-06, 2.877938264098797e-05, -3.8888124247418959e-06, -9.182484129984373e-08, 3.9073235571712722e-06, -1.0590135696371576e-05, -1.1814031371030429e-06
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
    993.74083954914227, 1.4552699836529146, -0.48971417606855078, 1.0057863507275291, 0.63905394851042974, 0.2822354346134604, 0.70155064775999798
};
static constexpr int a1_coeffs_n_annual = 3;
static constexpr double a1_coeffs[] = {
    0.26787420061650224, -0.17021116394048272, 0.017225314872679982, -0.010814702937179154, -0.0064318773737189515, 0.00072775379263137247, 0.010668499058709322
};
static constexpr int b1_coeffs_n_annual = 3;
static constexpr double b1_coeffs[] = {
    0.35467370463775821, -0.21650600633903522, 0.045200957367788751, -0.035868189250888749, 0.06547589613464469, -0.022932809113196787, -0.020332459084580229
};
static constexpr int a2_coeffs_n_annual = 3;
static constexpr double a2_coeffs[] = {
    -0.163817728726215, 0.12187984424469134, 0.053082073444650164, 0.062621002651114824, 0.044596637798804704, -0.045391331046554816, -0.020318575411327616
};
static constexpr int b2_coeffs_n_annual = 3;
static constexpr double b2_coeffs[] = {
    -0.45077744175117623, -0.025811414760719117, -0.023315524709236281, 0.038106942473252282, -0.028267636983741715, 0.021178874454994136, 0.0099458438450665456
};
static constexpr int a3_coeffs_n_annual = 3;
static constexpr double a3_coeffs[] = {
    0.017559707497729678, 0.11498580301615106, -0.038918433578680078, -0.022015452748607889, -0.011892212569845798, -0.0088119777880463386, 0.0026228622236016087
};
static constexpr int b3_coeffs_n_annual = 3;
static constexpr double b3_coeffs[] = {
    0.048642662615614222, 0.10348964927211025, 0.0091797079348757294, 0.021996988313875223, 0.01623740814723543, -0.0022105184157002597, 0.0094387798825203775
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
    1.1203657802575655, -0.040534694118315298, -0.010205323265297555,
    -0.0026878127006143698, 0.9691590588942125, 0.048732733892152026,
    -0.0070347493248891215, 0.0042050602860355182, 1.0241052308929819,
    -0.11450372087966018, 0.015177455153612684, 0.0096197630142992041,
    0.0030445951417294812, -0.016053235739626887, -0.028115986651337355,
    0.013780133381573847, 0.012261518322953732, -0.10086415119209266,
    -0.011546702492222935, 0.021594387219192157, 0.0089796523719676488,
    -0.0011057391989794281, 0.010270446434170466, -0.00069967069280300019,
    -0.0052230501782687677, 0.0050238774472401854, 0.04257235015639161
};
static constexpr double L[9] = {
    0.056797786579639083, 0, 0,
    0.0040319648445876028, 0.21693728849886487, 0,
    -0.044445238165427663, 0.006167442532934839, 0.18622070132556445
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