rwildcard	=	$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

COMPILER	=	clang++

################################################################################
# Desktop editor                                                               #
################################################################################
EDITORNAME		:=	lambda
EDITORPATH		:=	./$(EDITORNAME)
EDITORINC		:=	-I./src
EDITORLIBS		:=	-lm			\
					-ldl		\
					-lSDL2
EDITORDEFS	:=	-DUSE_SDL

SOURCE_CPP		:=	$(call rwildcard,./src/,*.cpp)
ALL_SOURCE 		:=	$(SOURCE_CPP)
EDITOR_OBJS		+=	$(SOURCE_CPP:.cpp=.o)

DEBUG_DEFINES	:=	-D_DEBUG
RELEASE_DEFINES :=	-DNDEBUG

################################################################################
# Compilation options                                                          #
################################################################################
CFLAGS = -O0 -g -c -gdwarf-4 -W -Wall -Werror -mavx -std=c++20 $(DEFINES)

.cpp.o:
	$(COMPILER) $(CFLAGS) $(INCLUDES) $< -o $@

.PHONY: all clean lambda

################################################################################
# Build targets                                                                #
################################################################################
lambda: INCLUDES = $(EDITORINC)
lambda: DEFINES = $(EDITORDEFS) $(DEBUG_DEFINES)
lambda: $(EDITOR_OBJS) $(ALL_SOURCE)
	$(COMPILER) $(EDITOR_OBJS) $(EDITORINC) $(EDITORLIBS) -o $(EDITORPATH)

# Cleaning rule
clean:
	$(RM) $(EDITOR_OBJS)
	$(RM) $(EDITORPATH)
	$(RM) *~

all: clean lambda
