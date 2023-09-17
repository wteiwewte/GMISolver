#ifndef GMISOLVER_MATRIXTYPES_H
#define GMISOLVER_MATRIXTYPES_H

#include <deque>
#include <vector>

enum class MatrixRepresentationType : int8_t {
  NORMAL = 0,
  SPARSE,
};

template <typename T> struct IndexedValue {
  T _data;
  int _index;
};

template <typename T,
          template <typename, typename...> class UnderlyingContainerT =
              std::deque,
          typename... Args>
struct SparseVector {
  const T &operator[](const int index) const { return _normalVec[index]; }

  UnderlyingContainerT<IndexedValue<T>, Args...> _indexedValues;
  std::vector<T> _normalVec;
};

template <typename T> struct SparseMatrixRepresentation {
  void clear() {
    _rows.clear();
    _columns.clear();
  }

  std::vector<SparseVector<T>> _rows;
  std::vector<SparseVector<T>> _columns;
};

template <typename T> using Matrix = std::vector<std::vector<T>>;

template <typename T> struct MatrixRepresentation {
  void clear() {
    _rows.clear();
    _columns.clear();
  }

  std::vector<std::vector<T>> _rows;
  std::vector<std::vector<T>> _columns;
};

template <typename T> struct ElementaryMatrix {
  std::vector<T> _vec;
  T _pivotingTermInverse{0.0};
  int _pivotRowIdx{0};
};

template <typename T> struct SparseElementaryMatrix {
  SparseVector<T> _sparseVec;
  T _pivotingTermInverse{0.0};
  int _pivotRowIdx{0};
};

#endif // GMISOLVER_MATRIXTYPES_H
