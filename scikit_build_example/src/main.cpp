#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <AudioFile.h>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace plt = matplot;
namespace py = pybind11;

struct Plot {
    std::vector<double> time;
    std::vector<double> values;
};


Plot audio_sample(const std::string& filePath) {
    Plot data;
    AudioFile<double> audioFile;

    if (!audioFile.load(filePath)) {
        std::cout << "Nie uda³o siê wczytaæ pliku audio." << std::endl;
        return data;
    }

    auto& samples = audioFile.samples[0];
    double sampleRate = audioFile.getSampleRate();

    int newNumSamples = (2 * sampleRate);
    if (samples.size() > newNumSamples) {
        samples.resize(newNumSamples);
    }
    data.values = samples;
    data.time.resize(samples.size());
    for (int i = 0; i < samples.size(); ++i) {
        data.time[i] = i / sampleRate;
    }

    return data;
}

Plot gen_signal(std::string n, double frequency) {
    Plot data;

    double pi = std::acos(-1);
    int sampling = 1000;
    double amplitude = 1;
    double phasor = 0;
    double step = frequency / sampling;
    double value = 0;

    for (double t = 0.0; t < 1.0; t += 1.0 / sampling) {
        if (n == "sin") value = amplitude * std::sin(2 * pi * frequency * t);
        else if (n == "cos")value = amplitude * std::cos(2 * pi * frequency * t);
        else if (n == "sqr") {
            if (phasor < 0.5) {
                value = amplitude;
            }
            else {
                value = -amplitude;
            }
        }
        else if (n == "saw") {
            if (t == 0.0) {
                value = 0.0;
            }
            else if (t != 0.0) {
                value += step;
                if (value >= 1.0) {
                    value -= 1.0;
                }
            }
        }
        phasor += frequency / sampling;
        if (phasor >= 1) phasor = 0.0;
        data.time.push_back(t);
        data.values.push_back(value);
    }
    return data;
}

void show_plot(const Plot& data) {
    if ((data.time).empty() || (data.values).empty()) {
        std::cout << "Wprowadzone bledne dane";
        return;
    }
    plt::plot(data.time, data.values);
    plt::xlabel("Czas [s]");
    plt::ylabel("Wartosc");
    plt::title("Wykres Y(t)");
    plt::show();
}

Plot diff(const Plot& data) {
    Plot diff_data;
    std::vector<double> diff;
    if (data.time.empty() || data.values.empty()) {
        std::cerr << "Brak danych do obliczenia pochodnej." << std::endl;
        return diff_data;
    }
    double timestep = data.time[1] - data.time[0];
    double derivative = 0.0;

    for (int i = 0; i < data.values.size(); i++) {
        if (i == 0) {
            derivative = (data.values[i + 1] - data.values[i]) / timestep;
        }
        else if (i == data.values.size() - 1) {
            derivative = (data.values[i] - data.values[i - 1]) / timestep;
        }
        else {
            derivative = (data.values[i + 1] - data.values[i - 1]) / timestep;
        }
        diff_data.time.push_back(data.time[i]);
        diff.push_back(derivative);
    }

    double min = *std::min_element(diff.begin(), diff.end());
    double max = *std::max_element(diff.begin(), diff.end());

    for (int i = 0; i < diff.size(); i++) {
        double scaled = 2 * (diff[i] - min) / (max - min) - 1;
        diff_data.values.push_back(scaled);
    }


    return diff_data;
}


Plot dft(const Plot& data) {
    Plot dft_data;
    double pi = std::acos(-1);
    if (data.time.empty() || data.values.empty()) {
        std::cout << "Brak danych do obliczenia DFT." << std::endl;
        return dft_data;
    }

    int N = data.values.size();


    for (int n = 0; n < N; n++) {
        std::complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; k++) {
            double angle = 2 * pi * k * n / N;
            sum += data.values[k] * std::exp(std::complex<double>(0.0, angle));
        }
        dft_data.time.push_back(data.time[n]);
        dft_data.values.push_back(std::abs(sum));
    }

    return dft_data;
}


Plot idft(const Plot& dft_data) {
    Plot idft_data;
    double pi = std::acos(-1);
    if (dft_data.time.empty() || dft_data.values.empty()) {
        std::cerr << "Brak danych do obliczenia IDFT." << std::endl;
        return idft_data;
    }

    int N = dft_data.values.size();
    idft_data.time.resize(N);
    idft_data.values.resize(N);

    for (int n = 0; n < N; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; ++k) {
            double angle = 2 * pi * k * n / N;
            double real_part = dft_data.values[k] * std::cos(angle) - dft_data.values[k] * std::sin(angle);
            sum += std::complex<double>(real_part, 0.0);
        }
        idft_data.time[n] = n;
        idft_data.values[n] = sum.real() / N;
    }
    return idft_data;
}

PYBIND11_MODULE(_core, m) {

    py::class_<Plot>(m, "Plot")
        .def(py::init<>())
        .def_readwrite("time", &Plot::time)
        .def_readwrite("values", &Plot::values);

    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("gen_signal", &gen_signal, R"pbdoc(
        Generate signal:
       "sin" - sinus
       "cos" - cosinus
       "sqr" - square 
       "saw" - sawtooth 
    )pbdoc");

    m.def("show_plot", &show_plot, R"pbdoc(
        Shows plot
    )pbdoc");

    m.def("audio_sample", &audio_sample, R"pbdoc(
        Generate audio signal
    )pbdoc");

    m.def("diff", &diff, R"pbdoc(
        Creates derivative of signal
    )pbdoc");

    m.def("dft", &dft, R"pbdoc(
        Creates dft of signal
    )pbdoc");


    m.def("idft", &idft, R"pbdoc(
        Creates idft of signal
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
