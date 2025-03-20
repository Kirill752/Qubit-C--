# Компилятор и флаги
CXX = g++
MESH_CXXFLAGS = -std=c++23 -O3 -IMesh/include
MESH_LDFLAGS = -lgmsh

# Исходные файлы
OBJ_DIR = obj

MESH_SRC_DIR = Mesh/src
MESH_DIR = Mesh

MESH_SRCS = $(wildcard $(MESH_SRC_DIR)/*.cpp)
MESH_OBJS = $(patsubst $(MESH_SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(MESH_SRCS))

# Исполняемый файл
MESH_TARGET = $(MESH_DIR)/mesh

# Правила сборки
all: mesh

mesh: $(MESH_TARGET)

$(MESH_TARGET): $(MESH_OBJS)
	$(CXX) $(MESH_CXXFLAGS) -o $(MESH_TARGET) Mesh/main.cpp $(MESH_OBJS) $(MESH_LDFLAGS)

$(OBJ_DIR)/%.o: $(MESH_SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(MESH_CXXFLAGS) -c $< -o $@

# Создание директории obj, если она не существует
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)


CALC_CXXFLAGS = -std=c++23 -O3 -I/usr/include/eigen3
CALC_LDFLAGS = -lmfem

CALC_DIR = Calculate

# Исполняемый файл
CALC_TARGET = $(CALC_DIR)/calculate
CALC_OBJS = $(OBJ_DIR)/Solver.o

calculate: $(CALC_TARGET)

$(CALC_TARGET): $(CALC_OBJS) $(CALC_DIR)/main.cpp
	$(CXX) $(CALC_CXXFLAGS) -o $(CALC_TARGET) $(CALC_DIR)/main.cpp $(CALC_OBJS) $(CALC_LDFLAGS)

$(OBJ_DIR)/%.o: $(CALC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CALC_CXXFLAGS) -c $< -o $@

# Очистка
clean:
	rm -f $(OBJ_DIR)/*.o

.PHONY: all mesh calculate clean