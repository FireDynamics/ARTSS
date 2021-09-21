#ifndef ARTSS_SOURCE_UNIFORMRANDOM_H_
#define ARTSS_SOURCE_UNIFORMRANDOM_H_


#include <random>
#include <iostream>

#include "../interfaces/IRandomField.h"


class UniformRandom: public IRandomField {
 public:
    UniformRandom(real range, real step_size) :
        m_seed(std::random_device()()),
        m_steps(static_cast<int>(range / step_size)),
        m_step_size(step_size) {
        init_dist();
    }

    UniformRandom(real range, real step_size, int seed) :
        m_seed(seed),
        m_steps(static_cast<int>(range / step_size)),
        m_step_size(step_size) {
        init_dist();
    }

    ~UniformRandom() override = default;

    Field random_field(size_t size) override;

 private:
     void init_dist() {
         m_mt = std::mt19937(m_seed);
         m_dist = std::uniform_int_distribution<int>(-m_steps, m_steps);
     }

     const int m_seed;
     const int m_steps;
     const real m_step_size;
     std::mt19937 m_mt;
     std::uniform_int_distribution<int> m_dist;
};

#endif /* ARTSS_SOURCE_UNIFORMRANDOM_H_ */
