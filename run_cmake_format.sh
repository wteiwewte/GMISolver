find .  -not -path "./cmake-build-debug/*" -not -path "./cmake-build-release/*" -iname "CMakeLists.txt" | xargs cmake-format -i
