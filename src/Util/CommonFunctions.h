#ifndef GMISOLVER_COMMONFUNCTIONS_H
#define GMISOLVER_COMMONFUNCTIONS_H

const auto isInteger = [](const auto val) {
  double integerPart;
  return std::modf(val, &integerPart) == 0.0;
};

#endif // GMISOLVER_COMMONFUNCTIONS_H
