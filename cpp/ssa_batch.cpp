#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
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

struct Endpoint {
    double extent = 0.0;
    double fni = 0.0;
    double q = 0.0, a = 0.0, p = 0.0, n = 0.0, m = 0.0, g = 0.0;
    bool failed = false;
};

static double ramp(double multiplier, double t, double tau) {
    if (multiplier <= 0.0) return 1.0;
    return std::exp(std::log(multiplier) * (1.0 - std::exp(-t / std::max(1e-6, tau))));
}

static Endpoint simulate(const Row& r, bool treated, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    long long Q = r.i("Q0"), A = r.i("A0"), P = r.i("P0"), N = r.i("N0"), M = r.i("M0"), G = r.i("G0");
    double F = static_cast<double>(G) * r.d("eff_mean");
    const double t_end = r.d("t_end"), tau = r.d("tau");
    const std::string model = r.s("model");
    double t = 0.0;
    bool failed = false;
    for (long long event = 0; event < 3000000 && t < t_end; ++event) {
        double ma = 1.0, mp = 1.0, mm = 1.0, mi = 1.0, me = 1.0, ms = 1.0;
        if (treated && model != "null") {
            if (model == "full" || model == "proliferation") {
                ma = ramp(r.d("tx_activation"), t, tau);
                mp = ramp(r.d("tx_prolif"), t, tau);
            }
            if (model == "full" || model == "maturation") {
                mm = ramp(r.d("tx_maturation"), t, tau);
                ms = ramp(r.d("tx_survival"), t, tau);
            }
            if (model == "full" || model == "integration") {
                mi = ramp(r.d("tx_integration"), t, tau);
                me = ramp(r.d("tx_eff"), t, tau);
            }
        }

        std::vector<double> a = {
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
        for (double x : a) total += std::max(0.0, x);
        if (!(total > 0.0) || !std::isfinite(total)) break;
        double u1 = std::max(unif(rng), std::numeric_limits<double>::min());
        t += -std::log(u1) / total;
        if (t > t_end) break;
        double target = unif(rng) * total, cumulative = 0.0;
        size_t reaction = a.size() - 1;
        for (size_t j = 0; j < a.size(); ++j) {
            cumulative += a[j];
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
                    --M; ++G;
                    double mean = std::max(1e-6, r.d("eff_mean") * me);
                    double cv = std::max(1e-6, r.d("eff_cv"));
                    double sigma2 = std::log(1.0 + cv * cv);
                    double mu = std::log(mean) - 0.5 * sigma2;
                    std::lognormal_distribution<double> ln(mu, std::sqrt(sigma2));
                    F += ln(rng);
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
    e.q = Q; e.a = A; e.p = P; e.n = N; e.m = M; e.g = G;
    e.failed = failed;
    return e;
}

static double mean(const std::vector<double>& x) {
    double s = 0.0; for (double v : x) s += v; return x.empty() ? std::nan("") : s / x.size();
}
static double sd(const std::vector<double>& x) {
    if (x.size() < 2) return 0.0;
    double m = mean(x), s = 0.0; for (double v : x) s += (v-m)*(v-m);
    return std::sqrt(s / (x.size()-1));
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: ssa_batch design.csv output.csv\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) throw std::runtime_error("cannot open design file");
    std::ofstream out(argv[2]);
    if (!out) throw std::runtime_error("cannot open output file");
    std::string line;
    std::getline(in, line);
    auto header = split_csv(line);
    out << "set_id,design_group,model,n_rep,failed_reps,control_extent,treat_extent,delta_extent,sd_delta_extent,control_fni,treat_fni,delta_fni,sd_delta_fni,treat_Q,treat_A,treat_P,treat_N,treat_M,treat_G\n";
    out << std::setprecision(10);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto cells = split_csv(line);
        if (cells.size() != header.size()) throw std::runtime_error("malformed design row");
        Row r; for (size_t i = 0; i < header.size(); ++i) r.v[header[i]] = cells[i];
        int reps = static_cast<int>(r.i("reps"));
        uint64_t base_seed = static_cast<uint64_t>(r.i("seed"));
        std::vector<double> ce, te, de, cf, tf, df, tq, ta, tp, tn, tm, tg;
        int failures = 0;
        for (int rep = 0; rep < reps; ++rep) {
            uint64_t s = mix64(base_seed ^ mix64(static_cast<uint64_t>(rep + 1)));
            Endpoint c = simulate(r, false, s);
            Endpoint t = simulate(r, true, s);
            failures += static_cast<int>(c.failed) + static_cast<int>(t.failed);
            ce.push_back(c.extent); te.push_back(t.extent); de.push_back(t.extent-c.extent);
            cf.push_back(c.fni); tf.push_back(t.fni); df.push_back(t.fni-c.fni);
            tq.push_back(t.q); ta.push_back(t.a); tp.push_back(t.p); tn.push_back(t.n); tm.push_back(t.m); tg.push_back(t.g);
        }
        out << r.s("set_id") << ',' << r.s("design_group") << ',' << r.s("model") << ',' << reps << ',' << failures << ','
            << mean(ce) << ',' << mean(te) << ',' << mean(de) << ',' << sd(de) << ','
            << mean(cf) << ',' << mean(tf) << ',' << mean(df) << ',' << sd(df) << ','
            << mean(tq) << ',' << mean(ta) << ',' << mean(tp) << ',' << mean(tn) << ',' << mean(tm) << ',' << mean(tg) << '\n';
    }
    return 0;
}
