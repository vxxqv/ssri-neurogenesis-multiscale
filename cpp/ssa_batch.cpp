#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct Row {
    std::unordered_map<std::string, std::string> v;
    double d(const std::string& key) const { return std::stod(v.at(key)); }
    long long i(const std::string& key) const { return std::stoll(v.at(key)); }
    std::string s(const std::string& key) const { return v.at(key); }
};

struct Endpoint {
    double extent = 0.0;
    double fni = 0.0;
    double q = 0.0;
    double a = 0.0;
    double p = 0.0;
    double n = 0.0;
    double m = 0.0;
    double g = 0.0;
    bool failed = false;
};

static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::string cell;
    bool quote = false;
    for (char c : line) {
        if (c == '"') quote = !quote;
        else if (c == ',' && !quote) { out.push_back(cell); cell.clear(); }
        else cell.push_back(c);
    }
    out.push_back(cell);
    return out;
}

static uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static double ramp(double multiplier, double t, double tau) {
    return std::exp(std::log(std::max(1e-9, multiplier)) * (1.0 - std::exp(-t / std::max(1e-6, tau))));
}

static Endpoint simulate(const Row& r, bool treated, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    long long Q = r.i("Q0"), A = r.i("A0"), P = r.i("P0"), N = r.i("N0"), M = r.i("M0"), G = r.i("G0");
    double F = static_cast<double>(G) * r.d("eff_mean");
    const double t_end = r.d("t_end"), tau = r.d("tau");
    double t = 0.0;
    bool failed = false;
    for (long long event = 0; event < 3000000 && t < t_end; ++event) {
        double ma = 1.0, mp = 1.0, mm = 1.0, mi = 1.0, me = 1.0, ms = 1.0;
        if (treated) {
            if (r.i("activation_on")) ma = ramp(r.d("tx_activation"), t, tau);
            if (r.i("proliferation_on")) mp = ramp(r.d("tx_prolif"), t, tau);
            if (r.i("maturation_on")) mm = ramp(r.d("tx_maturation"), t, tau);
            if (r.i("integration_on")) mi = ramp(r.d("tx_integration"), t, tau);
            if (r.i("efficacy_on")) me = ramp(r.d("tx_eff"), t, tau);
            if (r.i("survival_on")) ms = ramp(r.d("tx_survival"), t, tau);
        }
        std::vector<double> propensities = {
            r.d("k_qa") * ma * Q,
            r.d("k_aq") * A,
            r.d("b_a") * mp * A,
            r.d("d_a") * A,
            r.d("b_p") * mp * P,
            r.d("k_pn") * P,
            r.d("d_p") * P,
            r.d("k_nm") * mm * N,
            r.d("d_n") * ms * N,
            r.d("k_mg") * mi * M,
            r.d("d_m") * ms * M,
            r.d("d_g") * ms * G
        };
        double total = 0.0;
        for (double value : propensities) total += std::max(0.0, value);
        if (!(total > 0.0) || !std::isfinite(total)) break;
        double u1 = std::max(unif(rng), std::numeric_limits<double>::min());
        t += -std::log(u1) / total;
        if (t > t_end) break;
        double target = unif(rng) * total, cumulative = 0.0;
        size_t reaction = propensities.size() - 1;
        for (size_t j = 0; j < propensities.size(); ++j) {
            cumulative += propensities[j];
            if (target <= cumulative) { reaction = j; break; }
        }
        switch (reaction) {
            case 0: if (Q > 0) { --Q; ++A; } break;
            case 1: if (A > 0) { --A; ++Q; } break;
            case 2: if (A > 0) ++P; break;
            case 3: if (A > 0) --A; break;
            case 4: if (P > 0) ++P; break;
            case 5: if (P > 0) { --P; ++N; } break;
            case 6: if (P > 0) --P; break;
            case 7: if (N > 0) { --N; ++M; } break;
            case 8: if (N > 0) --N; break;
            case 9:
                if (M > 0) {
                    --M;
                    ++G;
                    double expected = std::max(1e-6, r.d("eff_mean") * me);
                    double cv = std::max(1e-6, r.d("eff_cv"));
                    double sigma2 = std::log(1.0 + cv * cv);
                    double mu = std::log(expected) - 0.5 * sigma2;
                    std::lognormal_distribution<double> draw(mu, std::sqrt(sigma2));
                    F += draw(rng);
                }
                break;
            case 10: if (M > 0) --M; break;
            case 11:
                if (G > 0) { F -= F / static_cast<double>(G); --G; F = std::max(0.0, F); }
                break;
        }
        long long total_cells = Q + A + P + N + M + G;
        if (total_cells < 0 || total_cells > 1000000) { failed = true; break; }
    }
    Endpoint e;
    e.extent = static_cast<double>(M + G);
    e.fni = F;
    e.q = Q;
    e.a = A;
    e.p = P;
    e.n = N;
    e.m = M;
    e.g = G;
    e.failed = failed;
    return e;
}

static double mean(const std::vector<double>& values) {
    double total = 0.0;
    for (double value : values) total += value;
    return values.empty() ? std::nan("") : total / values.size();
}

static double sd(const std::vector<double>& values) {
    if (values.size() < 2) return 0.0;
    double center = mean(values), total = 0.0;
    for (double value : values) total += (value - center) * (value - center);
    return std::sqrt(total / (values.size() - 1));
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: ssa_batch design.csv aggregate.csv replicates.csv\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    std::ofstream aggregate(argv[2]);
    std::ofstream replicate(argv[3]);
    if (!in || !aggregate || !replicate) throw std::runtime_error("cannot open input or output file");
    std::string line;
    std::getline(in, line);
    auto header = split_csv(line);
    aggregate << "set_id,base_id,design_group,model,t_end,n_rep,failed_reps,control_extent,treat_extent,delta_extent,sd_delta_extent,control_fni,treat_fni,delta_fni,sd_delta_fni,treat_Q,treat_A,treat_P,treat_N,treat_M,treat_G\n";
    replicate << "set_id,base_id,design_group,model,t_end,replicate,failed,control_extent,treat_extent,delta_extent,control_fni,treat_fni,delta_fni,control_Q,control_A,control_P,control_N,control_M,control_G,treat_Q,treat_A,treat_P,treat_N,treat_M,treat_G\n";
    aggregate << std::setprecision(10);
    replicate << std::setprecision(10);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto cells = split_csv(line);
        if (cells.size() != header.size()) throw std::runtime_error("malformed design row");
        Row r;
        for (size_t i = 0; i < header.size(); ++i) r.v[header[i]] = cells[i];
        int reps = static_cast<int>(r.i("reps"));
        if (reps <= 0) continue;
        uint64_t base_seed = static_cast<uint64_t>(r.i("seed"));
        std::vector<double> ce, te, de, cf, tf, df, tq, ta, tp, tn, tm, tg;
        int failures = 0;
        for (int rep = 0; rep < reps; ++rep) {
            uint64_t seed = mix64(base_seed ^ mix64(static_cast<uint64_t>(rep + 1)));
            Endpoint control = simulate(r, false, seed);
            Endpoint treatment = simulate(r, true, seed);
            bool failed = control.failed || treatment.failed;
            failures += static_cast<int>(failed);
            ce.push_back(control.extent);
            te.push_back(treatment.extent);
            de.push_back(treatment.extent - control.extent);
            cf.push_back(control.fni);
            tf.push_back(treatment.fni);
            df.push_back(treatment.fni - control.fni);
            tq.push_back(treatment.q);
            ta.push_back(treatment.a);
            tp.push_back(treatment.p);
            tn.push_back(treatment.n);
            tm.push_back(treatment.m);
            tg.push_back(treatment.g);
            replicate << r.s("set_id") << ',' << r.s("base_id") << ',' << r.s("design_group") << ',' << r.s("model") << ',' << r.d("t_end") << ',' << rep + 1 << ',' << failed << ','
                      << control.extent << ',' << treatment.extent << ',' << treatment.extent - control.extent << ','
                      << control.fni << ',' << treatment.fni << ',' << treatment.fni - control.fni << ','
                      << control.q << ',' << control.a << ',' << control.p << ',' << control.n << ',' << control.m << ',' << control.g << ','
                      << treatment.q << ',' << treatment.a << ',' << treatment.p << ',' << treatment.n << ',' << treatment.m << ',' << treatment.g << '\n';
        }
        aggregate << r.s("set_id") << ',' << r.s("base_id") << ',' << r.s("design_group") << ',' << r.s("model") << ',' << r.d("t_end") << ',' << reps << ',' << failures << ','
                  << mean(ce) << ',' << mean(te) << ',' << mean(de) << ',' << sd(de) << ','
                  << mean(cf) << ',' << mean(tf) << ',' << mean(df) << ',' << sd(df) << ','
                  << mean(tq) << ',' << mean(ta) << ',' << mean(tp) << ',' << mean(tn) << ',' << mean(tm) << ',' << mean(tg) << '\n';
    }
    return 0;
}
