#ifndef ARTSS_INTERFACES_IRANDOM_H_
#define ARTSS_INTERFACES_IRANDOM_H_

#include "../field/Field.h"

class IRandomField {
 public:
     virtual ~IRandomField() = default;
     virtual Field random_field() = 0;
};


#endif
