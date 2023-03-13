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
CFLAGS = -O0 -g -gdwarf-4 -W -Wall -Werror -mavx $(DEBUG_DEFINES) $(DEFINES)
#CFLAGS = -O2 -W -Wall -Werror -mavx $(RELEASE_DEFINES) $(DEFINES)

.cpp.o:
	$(COMPILER) $(CFLAGS) -std=c++20 -c -stdlib=libc++ -fimplicit-modules -fimplicit-module-maps $(INCLUDES) $< -o $@

.PHONY: all clean lambda

################################################################################
# Build targets                                                                #
################################################################################
lambda: INCLUDES = $(EDITORINC)
lambda: DEFINES = $(EDITORDEFS)
lambda: $(EDITOR_OBJS) $(ALL_SOURCE)
	$(COMPILER) $(EDITOR_OBJS) $(EDITORINC) $(EDITORLIBS) -stdlib=libc++ -fimplicit-modules -fimplicit-module-maps -o $(EDITORPATH)

# Cleaning rule
clean:
	$(RM) $(EDITOR_OBJS)
	$(RM) $(EDITORPATH)
	$(RM) *~

all: clean lambda
