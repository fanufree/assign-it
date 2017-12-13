/*
 * XplorFormat.cpp
 *
 *  Created on: 2013-08-06
 *      Author: e4k2
 */

#include "XplorFormat.h"
#include "Residue.h"

tr1::unordered_map< int, tr1::unordered_map<string,string> > XplorFormat::pdb2psf;
bool XplorFormat::initFlag = false;

void XplorFormat::init()
{
	if (initFlag)
		return;

	tr1::unordered_map<string,string> ala;
	tr1::unordered_map<string,string> arg;
	tr1::unordered_map<string,string> asp;
	tr1::unordered_map<string,string> asn;
	tr1::unordered_map<string,string> cys;
	tr1::unordered_map<string,string> glu;
	tr1::unordered_map<string,string> gln;
	tr1::unordered_map<string,string> gly;
	tr1::unordered_map<string,string> his;
	tr1::unordered_map<string,string> ile;
	tr1::unordered_map<string,string> leu;
	tr1::unordered_map<string,string> lys;
	tr1::unordered_map<string,string> met;
	tr1::unordered_map<string,string> phe;
	tr1::unordered_map<string,string> pro;
	tr1::unordered_map<string,string> ser;
	tr1::unordered_map<string,string> thr;
	tr1::unordered_map<string,string> trp;
	tr1::unordered_map<string,string> tyr;
	tr1::unordered_map<string,string> val;

	arg["HB3"] = "HB1";
	arg["HG3"] = "HG1";
	arg["HD3"] = "HD1";
	asp["HB3"] = "HB1";
	asn["HB3"] = "HB1";
	cys["HB3"] = "HB1";
	glu["HB3"] = "HB1";
	glu["HG3"] = "HG1";
	gln["HB3"] = "HB1";
	gln["HG3"] = "HG1";
	gly["HA3"] = "HA1";
	his["HB3"] = "HB1";
	ile["HG13"] = "HG11";
	leu["HB3"] = "HB1";
	lys["HB3"] = "HB1";
	lys["HG3"] = "HG1";
	lys["HD3"] = "HD1";
	lys["HE3"] = "HE1";
	met["HB3"] = "HB1";
	met["HG3"] = "HG1";
	phe["HB3"] = "HB1";
	pro["H2"] = "HT2";
	pro["H3"] = "HT1";
	pro["HB3"] = "HB1";
	pro["HG3"] = "HG1";
	pro["HD3"] = "HD1";
	ser["HB3"] = "HB1";
	trp["HB3"] = "HB1";
	tyr["HB3"] = "HB1";

	pdb2psf[Residue::A] = ala;
	pdb2psf[Residue::R] = arg;
	pdb2psf[Residue::D] = asp;
	pdb2psf[Residue::N] = asn;
	pdb2psf[Residue::C] = cys;
	pdb2psf[Residue::E] = glu;
	pdb2psf[Residue::Q] = gln;
	pdb2psf[Residue::G] = gly;
	pdb2psf[Residue::H] = his;
	pdb2psf[Residue::I] = ile;
	pdb2psf[Residue::L] = leu;
	pdb2psf[Residue::K] = lys;
	pdb2psf[Residue::M] = met;
	pdb2psf[Residue::F] = phe;
	pdb2psf[Residue::P] = pro;
	pdb2psf[Residue::S] = ser;
	pdb2psf[Residue::T] = thr;
	pdb2psf[Residue::W] = trp;
	pdb2psf[Residue::Y] = tyr;
	pdb2psf[Residue::V] = val;

	for (tr1::unordered_map< int, tr1::unordered_map<string,string> >::iterator it = pdb2psf.begin(); it != pdb2psf.end(); ++it)
	{
		int type = it->first;
		tr1::unordered_map<string,string>& mappings = it->second;
		if (type != Residue::PRO)
			mappings["H"] = "HN";
		mappings["H1"] = "HT1";
		mappings["H2"] = "HT2";
		mappings["H3"] = "HT3";
		mappings["OXT"] = "OT2";
	}

	initFlag = true;
}

string XplorFormat::translateProton(Residue* res, int resNumOffset, bool isMethyl, Atom* hx)
{
	if (!initFlag)
		init();

	tr1::unordered_map< int, tr1::unordered_map<string,string> >::iterator toPSF1 = XplorFormat::pdb2psf.find(res->type);
	tr1::unordered_map<string,string>& resMap = toPSF1->second;
	string name = hx->name;

	// translate amide H of residue 1 to HT# if hx->name == H
	if (res->num+resNumOffset == 1 && name == "H")
		return "HT#";

	tr1::unordered_map<string,string>::iterator nameMap = resMap.find(name);
	if (nameMap != resMap.end())
	{
		string temp = nameMap->second;
		if (isMethyl)
			return temp.substr(0,temp.length()-1)+'#';
		else
			return temp;
	}
	else
	{
		if (isMethyl)
			return name.substr(0,name.length()-1)+'#';
		else
			return name; // else name stays the same
	}
	return name;
}
