CXX=g++
LD=g++

CXXFLAGS += -Wall -std=c++11 -g $(shell root-config --cflags) 
CXXFLAGS += -I$(WCSIMDIR)/include
CXXFLAGS += -I./inc
LDFLAGS += $(shell root-config --ldflags) $(shell root-config --libs) -lTreePlayer -lMinuit2 
LDFLAGS += -L$(WCSIMDIR) -l WCSimRoot

.PHONY: clean Execs

OBJS1=MCNeuCapManager.o\
	 HitsManager.o HitCluster.o\
	 VertexFit.o ResTFcn.o\
	 NtagData.o NtagUtil.o

SRCDIR = ./src
OBJDIR = ./obj

SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJS = $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cc=.o)))

ntag: $(OBJS)
	$(RM) $@
	$(LD) $^ $(LDFLAGS) -o $@
	$(RM) $^
	@echo "$@ done"

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@[ -d $(OBJDIR) ]
	$(LD) $(CXXFLAGS) -o $@ -c $<

clean:
	$(RM) $(OBJDIR)/*.o
