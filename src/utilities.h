/**
 * @file utilities.h
 * @brief Utilties for Edwards Curves and points on them
 * @author Graham Enos
 *
 * This file contains various utilities for e2c2.
 */


#ifndef _UTILITIES_H
#define _UTILITIES_H


#include <algorithm>        // For reverse
#include <sstream>          // For stringstream set_parameter hackery
#include <stdexcept>        // Exceptions


/// Namespace for our library
namespace e2c2 {

    /**
     * @p set_parameter(T& param, const std::string& value,
     *      const bool& hex_and_rev=false)
     * @brief Generic hackery to set a parameter to a string
     *
     * @tparam T        type of parameter
     *
     * This function sets the parameter "param" to the value given in the
     * string "value." Basically this is to smooth over some of the rough edges
     * of NTL.
     */
    template <class T>
    void set_parameter(T& param, const std::string& value,
            const bool& hex_and_rev=false) {
        std::stringstream ss;
        std::ostream& out = ss;
        std::istream& in = ss;
        std::string v(value.begin(), value.end());
        if (hex_and_rev) {
            reverse(v.begin(), v.end());
            out << "0x";
        }
        out << v;
        in >> param;
    }


    /**
     * @p InvalidParametersException
     * @brief Custom exception to be thrown when building a curve with invalid
     * parameters
     */
    class InvalidParametersException : public std::invalid_argument {
    public:
        InvalidParametersException() :
            std::invalid_argument("INVALID PARAMETERS") {}
    };


    /**
     * @p NotImplementedException
     * @brief Custom exception to be thrown when we reach the limits of current
     * implementation
     */
    class NotImplementedException : public std::runtime_error {
    public:
        NotImplementedException() :
            std::runtime_error("NOT YET IMPLEMENTED") {}
    };


    /**
     * @p DifferentCurvesException
     * @brief Custom exception to be thrown when attempting to operate on
     * points from different curves
     */
    class DifferentCurvesException : public std::invalid_argument {
    public:
        DifferentCurvesException() :
            std::invalid_argument("THESE POINTS BELONG TO DIFFERENT CURVES") {}
    };
}
#endif  // _UTILITIES_H
