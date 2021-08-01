//
// Created by Prakhar Srivastav (https://github.com/prakhar1989/progress-cpp)
// Output slightly modified to include more brackets and iterations / second by Alec Glisman on
// 1-5-12

#ifndef PROGRESSBAR_PROGRESSBAR_HPP
#define PROGRESSBAR_PROGRESSBAR_HPP

// External dependencies
#include <chrono>   // Timing parts of the program
#include <cmath>    // math operations, such as ceil() and floor()
#include <iomanip>  // String formatting (std::setprecision())
#include <iostream> // data printing

class ProgressBar
{
  private:
    unsigned int ticks = 0;

    const unsigned int                          total_ticks;
    const unsigned int                          bar_width;
    const char                                  complete_char   = '=';
    const char                                  incomplete_char = ' ';
    const std::chrono::steady_clock::time_point start_time      = std::chrono::steady_clock::now();

  public:
    ProgressBar(unsigned int total, unsigned int width, char complete, char incomplete)
        : total_ticks{total}, bar_width{width}, complete_char{complete}, incomplete_char{incomplete}
    {
    }

    ProgressBar(unsigned int total, unsigned int width) : total_ticks{total}, bar_width{width}
    {
    }

    unsigned int
    operator++()
    {
        return ++ticks;
    }

    void
    display() const
    {
        float        progress = (float)ticks / total_ticks;
        unsigned int pos      = (int)(bar_width * progress);

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        auto                                  time_elapsed =
            std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();

        std::cout << "[";

        for (unsigned int i = 0; i < bar_width; ++i)
        {
            if (i < pos)
                std::cout << complete_char;
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << incomplete_char;
        }
        double hourElapsed = std::floor(double(time_elapsed) / 3600.0);
        double minElapsed  = static_cast<int>(double(time_elapsed) / 60.0) % 60;
        double secElapsed  = static_cast<int>(double(time_elapsed)) % 60;

        double iterationPerSecond = double(ticks) / double(time_elapsed);

        double hourLeft = std::floor(double(total_ticks - ticks) / (iterationPerSecond * 3600));
        double minLeft =
            static_cast<int>(double(total_ticks - ticks) / (iterationPerSecond * 60.0)) % 60;
        double secLeft = static_cast<int>(double(total_ticks - ticks) / iterationPerSecond) % 60;

        std::cout << "] " << int(progress * 100.0) << "% [" << std::setprecision(0) << std::fixed
                  << hourElapsed << ":" << std::setprecision(0) << std::fixed << minElapsed << ":"
                  << std::setprecision(0) << std::fixed << secElapsed << " <- "
                  << std::setprecision(0) << std::fixed << hourLeft << ":" << std::setprecision(0)
                  << std::fixed << minLeft << ":" << std::setprecision(0) << std::fixed << secLeft
                  << ", " << std::setprecision(0) << std::fixed << iterationPerSecond << "it/s"
                  << "]\r";
        std::cout.flush();
    }

    void
    done() const
    {
        display();
        std::cout << std::endl;
    }
};

#endif // PROGRESSBAR_PROGRESSBAR_HPP
