TARGET = main.out

HDRS = \
	   project/include

SRCS = \
	   project/src/matrix.cpp \
	   project/src/eigenvalue.cpp \
	   project/src/main.cpp

.PHONY: all clean

all: $(SRCS)
	$(CXX) -Wall -Wextra -Werror -I $(HDRS) -o $(TARGET) $(CXXFLAGS) $(SRCS) 
	./$(TARGET)
clean:
	rm -rf $(TARGET)