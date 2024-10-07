#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

std::vector<std::vector<double> > read_from_file(const std::string &name) {
    std::ifstream file(name);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file");
    }
    std::string s;
    size_t i = 0;
    while (std::getline(file, s) && s[0] != '#') {
        ++i;
    }

    std::vector<std::vector<double> > data;
    while (std::getline(file, s)) {
        std::istringstream ss(s);
        std::vector<double> row;
        double value;
        while (ss >> value) {
            row.push_back(value);
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    file.close();
    return data;
}

void write_to_file(const std::vector<std::vector<double> > &data, const std::string &name) {
    std::ofstream outFile(name);
    if (!outFile) {
        return;
    }
    for (const auto &elem: data) {
        outFile << elem[0] << ' ' << elem[1] << ' ' << elem[2] << ' ' << elem[3] << std::endl;
    }
    outFile.close();
}


std::vector<std::vector<double> > mid_point(const std::vector<std::vector<double> > &data, const int p, int coord) {
    std::vector<int> index(p);
    for (int i = 0; i < p; ++i) {
        index[i] = static_cast<int>(i * (data.size() - 1) / (p - 1));
    }
    std::vector<std::vector<double> > finale = {{data[0][0], data[0][coord]}};
    for (int edge = 0; edge < p - 1; ++edge) {
        const double big = (*std::max_element(data.begin() + index[edge], data.begin() + index[edge + 1],
                                              [coord](const std::vector<double> &a, const std::vector<double> &b) { return a[coord] < b[coord]; }))[coord];
        const double small = (*std::min_element(data.begin() + index[edge], data.begin() + index[edge + 1],
                                                [coord](const std::vector<double> &a, const std::vector<double> &b) { return a[coord] < b[coord]; }))[coord];
        double centre = (big + small) / 2;
        auto element = *std::min_element(data.begin() + index[edge], data.begin() + index[edge + 1], [centre, coord](const std::vector<double> &a, const std::vector<double> &b) {
            return std::abs(a[coord] - centre) < std::abs(b[coord] - centre);
        });
        finale.push_back({element[0], element[coord]});
    }
    finale.push_back({data.back()[0], data.back()[coord]});
    return finale;
}

std::vector<std::vector<double> > linerSpline(const std::vector<std::vector<double> > &coords) {
    std::vector<std::vector<double> > cof;
    for (size_t edge = 0; edge < coords.size() - 1; ++edge) {
        double k = (coords[edge + 1][1] - coords[edge][1]) /
                   (coords[edge + 1][0] - coords[edge][0]);
        double b = coords[edge][1] - k * coords[edge][0];
        cof.push_back({b, k, coords[edge][0], coords[edge + 1][0]});
    }
    return cof;
}


double applySpline(const double p, const std::vector<std::vector<double> > &func) {
    for (const auto &elem: func) {
        if (elem[2] <= p && p <= elem[3]) {
            double sum = 0;
            for (int i = 0; i < elem.size() - 2; ++i) {
                sum += elem[i] * pow(p, i);
            }
            return sum;
        }
    }
    throw std::runtime_error("Point is outside the interpolation range.");
}

std::vector<std::vector<double> > LR(const std::vector<std::vector<double> > &l, const std::vector<std::vector<double> > &r, const std::vector<std::vector<double> > &f_real,
                                     const std::vector<std::vector<double> > &f_imag) {
    std::vector<std::vector<double> > finale;
    for (int i = 0; i < l.size(); ++i) {
        double g_real = 1, g_imag = 1, s22_r = 0, s22_i = 0, s1221_r = 0, s1221_i = 0;
        for (int j = 0; j < 4; ++j) {
            double s = sqrt((g_imag * (r[i][1] - applySpline(l[i][0], f_real)) * (1 - s22_r * g_real)) / (g_real * (r[i][8] - applySpline(l[i][0], f_imag)) * (1 - s22_i * g_imag)));
            s22_r = s * (l[i][8] - applySpline(l[i][0], f_imag)) / l[i][6];
            s22_i = 1 / s * (l[i][1] - applySpline(l[i][0], f_real)) / l[i][3];
            s1221_r = l[i][3] * s * (1 - s22_r * s22_i);
            s1221_i = l[i][6] * (1 / s) * (1 - s22_r * s22_i);
            g_real = (r[i][1] - applySpline(l[i][0], f_real)) / (r[i][1] * s22_r - applySpline(l[i][0], f_real) * s22_r + s1221_r);
            g_imag = (r[i][8] - applySpline(l[i][0], f_imag)) / (r[i][8] * s22_i - applySpline(l[i][0], f_imag) * s22_i + s1221_i);
        }
        finale.push_back({s22_r, s22_i, s1221_r, s1221_i});
    }
    return finale;
}

std::vector<std::vector<double> > LRT(const std::vector<std::vector<double> > &l, const std::vector<std::vector<double> > &r, const std::vector<std::vector<double> > &t,
                                      const std::vector<std::vector<double> > &f_real, const std::vector<std::vector<double> > &f_imag) {
    std::vector<std::vector<double> > finale;
    for (int i = 0; i < l.size(); ++i) {
        double g_real = 1, g_imag = 1, s22_r = 0, s22_i = 0, s1221_r = 0, s1221_i = 0;
        for (int j = 0; j < 4; ++j) {
            double s = sqrt((r[i][1] - applySpline(l[i][0], f_real)) * (1 - s22_r * g_real) / ((r[i][8] - applySpline(l[i][0], f_imag)) * (1 - s22_i * g_imag)));
            s22_r = s * (t[i][8] - applySpline(l[i][0], f_imag)) / t[i][6];
            s22_i = 1 / s * (t[i][1] - applySpline(l[i][0], f_real)) / t[i][3];
            s1221_r = t[i][3] * s * (1 - s22_r * s22_i);
            s1221_i = t[i][6] * (1 / s) * (1 - s22_r * s22_i);
            g_real = (r[i][1] - applySpline(l[i][0], f_real)) / (r[i][1] * s22_r - applySpline(l[i][0], f_real) * s22_r + s1221_r);
            g_imag = (r[i][8] - applySpline(l[i][0], f_imag)) / (r[i][8] * s22_i - applySpline(l[i][0], f_imag) * s22_i + s1221_i);
        }
        finale.push_back({s22_r, s22_i, s1221_r, s1221_i});
    }
    return finale;
}

std::vector<std::vector<double> > express(const std::vector<std::vector<double> > &l, const std::vector<std::vector<double> > &f_real, const std::vector<std::vector<double> > &f_imag) {
    std::vector<std::vector<double> > finale;
    for (int i = 0; i < l.size(); ++i) {
        double s22_r = 0, s22_i = 0, s1221_r = 0, s1221_i = 0;
        double s = sqrt(1);
        s22_r = s * (l[i][8] - applySpline(l[i][0], f_imag)) / l[i][6];
        s22_i = 1 / s * (l[i][1] - applySpline(l[i][0], f_real)) / l[i][3];
        s1221_r = l[i][3] * s * (1 - s22_r * s22_i);
        s1221_i = l[i][6] * (1 / s) * (1 - s22_r * s22_i);
        finale.push_back({s22_r, s22_i, s1221_r, s1221_i});
    }
    return finale;
}

int main() {
    std::vector<std::vector<double> > L = read_from_file("L120.s2p");
    std::vector<std::vector<double> > R = read_from_file("R120.s2p");
    std::vector<std::vector<double> > T = read_from_file("T30.s2p");
    auto sa = mid_point(L, 8, 1);
    auto sb = mid_point(L, 9, 2);
    auto sa_func = linerSpline(sa);
    auto sb_func = linerSpline(sb);
    auto lrt = LRT(L, R, T, sa_func, sb_func);
    auto lr = LR(L, R, sa_func, sb_func);
    auto expr = express(L, sa_func, sb_func);
    write_to_file(lrt, "lrt.txt");
    write_to_file(lr, "lr.txt");
    write_to_file(expr, "express.txt");
}
