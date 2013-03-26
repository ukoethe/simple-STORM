macro(ADD_RESOURCE RESOURCE_FILE HEADER_FILE VARIABLE_NAME SRCS)
    set(SRC "${HEADER_FILE}_.c")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/rc.cpp
"#include <fstream>
#include <iostream>
#include <locale>

int main(int argc, char* argv[])
{
    if (argc < 5)
        return 1;
    std::ifstream in(argv[1]);
    std::ofstream header(argv[2]);
    std::ofstream cfile(argv[3]);
    if (!in.is_open() || header.bad() || cfile.bad())
        return 1;
    std::string includeguard(argv[4]);

    std::locale loc;
    for (char &c : includeguard) {
        c = std::toupper(c, loc);
    }
    includeguard.append(\"_H\");

    header <<
        \"/* This is an automatically generated file. All changes\" << std::endl <<
        \"   will be lost */\" << std::endl <<
        \"#ifndef \" << includeguard << std::endl <<
        \"#define \" << includeguard << std::endl <<
        \"extern char \" << argv[4] << \"[];\" << std::endl <<
        \"#endif\";
    header.close();

    char cin;
    cfile <<
        \"/* This is an automatically generated file. All changes\" << std::endl <<
        \"   will be lost */\" << std::endl <<
        \"#include \\\"\" << argv[2] << \"\\\"\" << std::endl <<
        \"char \" << argv[4] << \"[] = {\";
    while (!in.eof()) {
        in.get(cin);
        cfile << (int)cin << ',';
    }
    cfile << \"0};\" << std::endl;
    cfile.close();
    in.close();
    return 0;
}
")
    add_executable(rc ${CMAKE_CURRENT_BINARY_DIR}/rc.cpp)

    add_custom_command(OUTPUT ${SRC} COMMAND rc ARGS "${RESOURCE_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${HEADER_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${SRC}" ${VARIABLE_NAME} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS rc)
    set(${SRCS} ${${SRCS}} ${SRC})
endmacro(ADD_RESOURCE)
