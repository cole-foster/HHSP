
# main- build and search
add_executable(main main.cpp)
target_link_libraries(main src)
target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/include)
install(TARGETS main)

# optimization of GHSP 
if(GHSP_OPT)
    add_executable(opt main-optimization.cpp)
    target_link_libraries(opt PRIVATE src)
    target_include_directories(opt PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_include_directories(opt PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(opt PRIVATE ${NLOPT})
    target_include_directories(opt PRIVATE ${NLOPT_DIRECTORIES})
    install(TARGETS opt)
endif()

# main- test
add_executable(test main-test.cpp)
target_link_libraries(test src)
target_include_directories(test PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(test PRIVATE ${PROJECT_SOURCE_DIR}/include)
install(TARGETS test)


