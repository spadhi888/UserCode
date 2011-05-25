
#ifndef DATASOURCE_H
#define DATASOURCE_H

#include "TString.h"
#include "TColor.h"
#include <iostream>

typedef UInt_t sources_t;

class DataSource {

	public:
		DataSource();
		DataSource(TString name, sources_t source, Color_t color = 0);
		~DataSource();

		TString         getName()       { return sourceName_; }
		sources_t       getSource()     { return source_; }
		sources_t 	getBit()	{ return 1ll << source_; }
		Color_t		getColor()	{ return color_; }
	private:
		TString         sourceName_;
		sources_t       source_;
		Color_t		color_;

};

#endif

