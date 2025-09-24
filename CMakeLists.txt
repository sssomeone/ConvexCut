cmake_minimum_required(VERSION 4.0)
project(ConvexCut)

set(CMAKE_CXX_STANDARD 20)

# -------------------- 第三方库 --------------------
find_package(Boost REQUIRED)            # Boost库
find_package(CGAL REQUIRED)             # 提供 CGAL::CGAL
find_package(Eigen3 3.3 REQUIRED)       # 提供变量 EIGEN3_INCLUDE_DIR

# If your main file is in a different location, update the path:
# add_executable(ConvexCut src/main.cpp)
# or if you have multiple source files:
# add_executable(ConvexCut main.cpp other_file.cpp)

add_executable(ConvexCut main.cpp)

# 设置包含目录
target_include_directories(ConvexCut
        PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# 链接库
target_link_libraries(ConvexCut
        PRIVATE
        CGAL::CGAL
        Boost::boost
)