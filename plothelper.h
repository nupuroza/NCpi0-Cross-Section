#ifndef PLOTHELPER_H
#define PLOTHELPER_H

template <typename T>
std::string to_string_prec(const T a_value, const int n)
{
    std::ostringstream out;
    out <<std::fixed<< std::setprecision(n) << a_value;
    return out.str();
}

#endif
