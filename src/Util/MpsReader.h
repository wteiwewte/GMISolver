#ifndef GMISOLVER_MPSREADER_H
#define GMISOLVER_MPSREADER_H

#include "src/DataModel/EnumTypes.h"
#include "src/DataModel/LinearProgram.h"

#include <map>
#include <optional>
#include <string>
#include <vector>

template <typename T> struct MpsReader {
public:
  std::optional<LinearProgram<T>> read(const std::string &filePath);

private:
  constexpr static int OBJECTIVE_ROW_CONSTRAINT_IDX = 0;

  std::optional<int> setupObjectiveRowVar(const bool currentSectionIsInteger,
                                          LinearProgram<T> &linearProgram);
  bool handleRowsSection(const std::string &readLine,
                         const std::vector<std::string> &lineParts,
                         LinearProgram<T> &linearProgram);
  bool handleColumnsSection(const std::string &readLine,
                            const std::vector<std::string> &lineParts,
                            LinearProgram<T> &linearProgram,
                            bool &currentSectionIsInteger);
  bool handleRHSSection(const std::string &readLine,
                        const std::vector<std::string> &lineParts,
                        LinearProgram<T> &linearProgram);
  bool handleBoundsSection(const std::string &readLine,
                           const std::vector<std::string> &lineParts,
                           LinearProgram<T> &linearProgram);
  bool handleRangesSection(const std::string &readLine,
                           const std::vector<std::string> &lineParts,
                           LinearProgram<T> &linearProgram);

  bool addNewRowInfo(const std::string &rowLabelStr, const RowType rowType,
                     LinearProgram<T> &linearProgram);
  std::optional<int> tryAddNewVar(const std::string &variableLabelStr,
                                  const bool currentSectionIsInteger,
                                  LinearProgram<T> &linearProgram);
  bool updateLpMatrix(const std::string &rowLabelStr,
                      const std::string &coefficientValueStr,
                      const int variableIdx, LinearProgram<T> &linearProgram);

  static void setBoundsForNonFreeVars(LinearProgram<T> &linearProgram);
  static void roundBoundsForIntegerVars(LinearProgram<T> &linearProgram);

  std::map<std::string, std::optional<int>> _rowLabelToRowIdxMap;
  std::map<std::string, int> _variableLabelToVariableIdxMap;
};

#endif // GMISOLVER_MPSREADER_H
