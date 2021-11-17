
#include "Settings.h"

std::string Settings::get(std::string path) {
    if (m_proxy.find(path) == m_proxy.end()) {
        auto val = m_params->get(path);
#ifndef BENCHMARKING
        m_logger->debug("first time reading: \"{}\" with value: \"{}\"", path, val);
#endif
        set(path, val);
    }

    auto iter = m_proxy.find(path);
#ifndef BENCHMARKING
     m_logger->debug("reading: \"{}\" with value: \"{}\"", path, iter->second);
#endif
    return iter->second;
}

