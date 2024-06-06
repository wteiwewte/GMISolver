#ifndef GMISOLVER_OVERLOADED_H
#define GMISOLVER_OVERLOADED_H

template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};

template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

#endif // GMISOLVER_OVERLOADED_H
