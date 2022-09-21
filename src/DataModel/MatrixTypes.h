#ifndef GMISOLVER_MATRIXTYPES_H
#define GMISOLVER_MATRIXTYPES_H

#include <deque>
#include <vector>

template <typename T>
struct SparseVector
{
  struct Element {
    T _data;
    int _index;
  };

  std::deque<Element> _elements;
  // TODO - use std::deque
  std::vector<T> _normalVec;
};

template <typename T>
struct SparseMatrixRepresentation
{
  void clear()
  {
    _rows.clear();
    _columns.clear();
  }

  std::vector<SparseVector<T>> _rows;
  std::vector<SparseVector<T>> _columns;
};

template <typename T> using Matrix = std::vector<std::vector<T>>;

template <typename T>
struct MatrixRepresentation
{
  void clear()
  {
    _rows.clear();
    _columns.clear();
  }

  std::vector<std::vector<T>> _rows;
  std::vector<std::vector<T>> _columns;
};

template <typename T>
struct ElementaryMatrix {
  std::vector<T> _vec;
  T _pivotingTermInverse{0.0};
  int _pivotRowIdx{0};
};

template <typename T>
struct SparseElementaryMatrix {
  SparseVector<T> _sparseVec;
  T _pivotingTermInverse{0.0};
  int _pivotRowIdx{0};
};

#endif // GMISOLVER_MATRIXTYPES_H
