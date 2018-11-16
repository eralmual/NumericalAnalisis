//
// Created by allan on 11/11/18.
//

#ifndef PROYECTO3_ISOLATIONCONTAINER_H
#define PROYECTO3_ISOLATIONCONTAINER_H

#include <string>
#include <Exception.hpp>
#include <vector>

namespace container {
    /**
     * This class stores the different isolation conditions of the plate
     */
    class Isolation {
    public:
        static bool topIsolated;
        static bool botIsolated;
        static bool leftIsolated;
        static bool rightIsolated;

        /**
         * This method runs as an extra_parser for the program_options, detects if -i is passed and parses
         * the borders which the users want isolated
         * @param s string to parse
         * @return pair of strings with the params (empty in this method)
         */
        static std::pair<std::string, std::string> parse_isolation(const std::string &s) {

            if (s.find("-i") == 0) {
                std::string rest = s.substr(2);
                for (char &it : rest) {
                    if (it == 't') {
                        Isolation::topIsolated = true;
                    } else if (it == 'b') {
                        Isolation::botIsolated = true;
                    } else if (it == 'l') {
                        Isolation::leftIsolated = true;
                    } else if (it == 'r') {
                        Isolation::rightIsolated = true;
                    } else {
                        throw anpi::Exception("isolation option not recognized");
                    }
                }

            }
            return std::make_pair(std::string(), std::string());
        }
        /**
         * this method stores the boolean parameters in a vector
         * @return vector containing the isolation for top bottom left right (in that order)
         */
        static std::vector<bool> convertToIsolationVector() {
            return std::vector<bool>{Isolation::topIsolated, Isolation::botIsolated, Isolation::leftIsolated,
                                     Isolation::rightIsolated};
        }
    };
}

namespace container {
    bool Isolation::topIsolated = false;
    bool Isolation::botIsolated = false;
    bool Isolation::leftIsolated = false;
    bool Isolation::rightIsolated = false;
}

#endif //PROYECTO3_ISOLATIONCONTAINER_H
