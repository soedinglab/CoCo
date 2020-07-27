//
// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
//

#include "Info.h"

int Info::verboseLevel = Info::INFO;

void Info::setVerboseLevel (int i) {
  verboseLevel = i;
}