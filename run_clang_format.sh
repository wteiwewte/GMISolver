find . -type f  -name "*.h" -not -path "*/cmake-build-debug/*" -not -path "*/cmake-build-release/*" | xargs clang-format -i
find . -type f  -name "*.cpp" -not -path "*/cmake-build-debug/*" -not -path "*/cmake-build-release/*" | xargs clang-format -i
