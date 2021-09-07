#ifndef ARTSS_INTERFACES_IRANDOM_H_
#define ARTSS_INTERFACES_IRANDOM_H_

#include "../field/Field.h"

class IRandomField {
 public:
     virtual ~IRandomField() = default;
     virtual Field random_field(size_t) = 0;
     Field random_field(const Field &in) {
         return random_field(in.get_size());
     }
};


#endif
