#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>


using namespace std;

int main() {
    // Initialization
    int size = 300;
    int size2 = 300;
    vector<vector<double>> T(size, vector<double>(size2, 0.0)); // Temperature profile
    vector<vector<double>> p(size, vector<double>(size2, 0.0)); // Phase field
    vector<vector<double>> theta(size, vector<double>(size2, 0.0));
    vector<vector<double>> grad_px(size, vector<double>(size2, 0.0));
    vector<vector<double>> grad_py(size, vector<double>(size2, 0.0));
    vector<vector<double>> eps(size, vector<double>(size2, 0.0));
    vector<vector<double>> eps_der(size, vector<double>(size2, 0.0));
    vector<vector<double>> grad_eps2x(size, vector<double>(size2, 0.0));
    vector<vector<double>> grad_eps2y(size, vector<double>(size2, 0.0));

    // Parameter values
    double del_t = 0.0001;
    double del_x = 0.03;
    double del_y = 0.03;
    double K = 2; // Dimensionless Latent Heat    
    double tou = 0.0003;
    double T_eq = 1;
    double gamma = 10;
    double alpha = 0.9;
    double a = 0.01;
    int J = 4; // Mode of anisotropy
    double delta = 0.05;
    double eps_bar = 0.01;

    vector<vector<double>> T_new(size, vector<double>(size2, 0.0));
    vector<vector<double>> p_new = p;
    double time = 1.4; // Total time for simulation
    int time_steps = static_cast<int>(time / del_t);

    // Initializing nucleus
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size2; ++j) {
            if ((i-size/2) * (i/size/2) + (j - size2 / 2) * (j - size2 / 2) < 20) {
                p[i][j] = 1;
            }
        }
    }

    double t = 0;
    int l = 0;
    vector<vector<double>> test(size, vector<double>(size2, 0.0));

    // Simulation starts
    while (l < time_steps) {
        for (int i = 1; i < size - 1; ++i) {
            for (int j = 1; j < size2 - 1; ++j) {
                grad_px[i][j] = (p[i + 1][j + 1] - p[i - 1][j + 1]) / (2 * del_x);
                grad_py[i][j] = (p[i + 1][j + 1] - p[i + 1][j - 1]) / (2 * del_y);

                if (grad_px[i][j] == 0) {
                    if (grad_py[i][j] > 0) {
                        theta[i][j] = 0.5 * M_PI;
                    }
                    else if (grad_py[i][j] < 0) {
                        theta[i][j] = 1.5 * M_PI;
                    }
                }
                if (grad_px[i][j] > 0) {
                    if (grad_py[i][j] < 0) {
                        theta[i][j] = 2 * M_PI + atan(grad_py[i][j] / grad_px[i][j]);
                    }
                    else if (grad_py[i][j] > 0) {
                        theta[i][j] = atan(grad_py[i][j] / grad_px[i][j]);
                    }
                }
                if (grad_px[i][j] < 0) {
                    theta[i][j] = M_PI + atan(grad_py[i][j] / grad_px[i][j]);
                }

                eps[i][j] = eps_bar * (1 + delta * cos(J * theta[i][j]));
                eps_der[i][j] = -eps_bar * J * delta * sin(J * theta[i][j]);
            }
        }

        for (int i = 1; i < size - 1; ++i) {
            for (int j = 1; j < size2 - 1; ++j) {
                double p11 = (p[i + 1][j + 1] - p[i + 1][j - 1]) / 2 * del_y;
                double p12= (p[i-1][j + 1] - p[i-1][j - 1]) / 2 * del_y;
                double part1 = -(eps[i + 1][j] * eps_der[i + 1][j] * p11 - eps[i - 1][j] * eps_der[i - 1][j] * p12) / 2 * del_x;

                double p21 = (p[i + 1][j + 1] - p[i - 1][j - 1]) / 2 * del_x;
                double p22 = (p[i+1][j-1] - p[i - 1][j - 1]) / 2 * del_x;
                double part2 = (eps[i][j+1] * eps_der[i][j+1] * p21 - eps[i][j-1] * eps_der[i][j-1] * p22) / 2 * del_y;

                double part3 = ((pow(eps[i + 1][j], 2) - pow(eps[i - 1][j], 2)) * (p[i + 1][j] - p[i - 1][j])) / pow(2 * del_x, 2) + ((pow(eps[i][j+1], 2) - pow(eps[i][j-1], 2)) * (p[i][j+1] - p[i][j-1])) / pow(2 * del_y, 2);
                double part4 = pow(eps[i][j], 2) * ((p[i + 1][j] + p[i - 1][j] - 2 * p[i][j]) / pow(del_x, 2) + (p[i][j + 1] + p[i][j - 1] - 2 * p[i][j]) / pow(del_y, 2));

                double m = alpha / M_PI * atan(gamma * (T_eq - T[i][j]));
                
                p_new[i][j] = p[i][j] + (part1 + part2 + part3 + part4 + p[i][j] * (1 - p[i][j]) * (p[i][j] - 0.5 + m)) * del_t / tou + a * p[i][j] * (1 - p[i][j]) * (static_cast<double>(rand()) / RAND_MAX - 0.5);
                T_new[i][j] = T[i][j] + ((T[i + 1][j] + T[i - 1][j] - 2 * T[i][j]) / (del_x * del_x) + (T[i][j + 1] + T[i][j - 1] - 2 * T[i][j]) / (del_y * del_y)) * del_t + K * (p_new[i][j] - p[i][j]);
                
                


            }
        }

        // Giving no flux boundary conditions
        for (int j = 0; j < size2; ++j) {
            p_new[0][j] = p_new[1][j];
            p_new[size - 1][j] = p_new[size - 2][j];
        }
        for (int i = 0; i < size; ++i) {
            p_new[i][0] = p_new[i][1];
            p_new[i][size2 - 1] = p_new[i][size2 - 2];
        }
        for (int j = 0; j < size2; ++j) {
            T_new[0][j] = T_new[1][j];
            T_new[size - 1][j] = T_new[size - 2][j];
        }
        for (int i = 0; i < size; ++i) {
            T_new[i][0] = T_new[i][1];
            T_new[i][size2 - 1] = T_new[i][size2 - 2];
        }

        T = T_new;
        p = p_new;
        t += del_t;
        l += 1;
        cout << l << endl;
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size2; ++j) {
            cout << p_new[i][j];
        }
    }

    //Save the data in a .csv file
    std::ofstream output_csv("output.csv");
    if (output_csv.is_open()) {
        for (const auto& row : p_new) {
            for (size_t i = 0; i < row.size(); ++i) {
                output_csv << row[i];
                if (i < row.size() - 1)
                    output_csv << ","; // Add a comma between values
            }
            output_csv << std::endl; // Newline for the next row
        }
        output_csv.close();
        std::cout << "CSV file created successfully." << std::endl;
    }
    else {
        std::cout << "Unable to create output file." << std::endl;
        return 1;
    }


    return 0;
}


