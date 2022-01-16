coil: coil.cpp
	g++ coil.cpp \
		-O2 \
		-I/usr/include/cairo -lcairo \
		`GraphicsMagick++-config --cppflags --cxxflags --ldflags --libs` \
		${MAGICKFLAGS} -o coil
