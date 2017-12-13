/*
 * PseudoContact.h
 *
 *  Created on: 2012-05-03
 *      Author: e4k2
 */

#ifndef PSEUDOCONTACT_H_
#define PSEUDOCONTACT_H_

#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include "Contact.h"

using namespace std;

/**
 *
 * Res1 Res2 {set of contacts between Res1 Res2}
 *
 * Atoms in the same residue with similar chemical shifts (according to some chemical shift tolerance)
 * are treated as a group (atom group). Contacts from such atoms (in residue A)
 * to atoms/atom group in another residue B are treated as a group of contacts known
 * as a PseudoContact.
 *
 * The order of the residues matter because NOEs have a direction
 */
class PseudoContact
{
public:
	tr1::unordered_set<Contact> contacts;
	int res1; // contact is from res1 to res2; if contacts is empty, then res1, res2 are < 0
	int res2;

	PseudoContact();
	PseudoContact(const Contact& c);
	PseudoContact(const PseudoContact& pc);
	// PseudoContact(PseudoContact&& pc);
	friend void swap(PseudoContact& first, PseudoContact& second);
	PseudoContact& operator=(PseudoContact pc);
	void add(const Contact& c); // no explicit check is made to ensure c is between the same pair of residues as those in contacts
	                        // the contact is not added if it already exists
	int numContacts() const;
	void print() const;
	virtual ~PseudoContact();
	bool operator==(const PseudoContact& pc) const; // contacts.size() must match as well as the contacts
	bool areSymmetric(const PseudoContact& pc) const; // true if contacts equals the reverse of pc.contacts; contacts.size() must match
	bool operator!=(const PseudoContact& pc) const;
	tr1::unordered_set<Contact>& getContacts();
	Contact getOneContact(); // returns one contact
	void getResPairOrdered(int& res1, int& res2) const; // returns the pair of residues that this PseudoContact represents; res1 <= res2
	                                                  // used for residue-level contacts
	void getResPair(int& res1, int& res2) const; // unlike getResPairOrdered, the ordering is not possibly changed; used for atomic level contacts
	                                       // where the ordering of the atom names must match the ordering of the residues
	void setReverse(); // reverses the direction of the contacts in contacts
	bool overlaps(const PseudoContact& pc) const; // returns true if contacts overlaps with at least one contact in pc
	                                     // it is possible that res1=res2' and res2=res1', in which case the reverse direction of the
	                                    // pc.contacts is used for comparison
	void addAll(const PseudoContact& pc); // adds all contacts in pc if not already in contacts; does not add anything if
	        // (res1 != pc.res1 || res2 != pc.res2) && (res1 != pc.res2 || res2 != pc.res1)
	int getSeqSep() const; // returns abs(res1-res2)
};

namespace std
{
	namespace tr1
	{
	template <>
		class hash<PseudoContact>
		{
			public:
			size_t operator()(const PseudoContact& c) const
			{
				unordered_set<Contact>::const_iterator it = c.contacts.begin();
				size_t ret = hash<Contact>()(*it);
				++it;
				for (; it != c.contacts.end(); ++it)
				{
					ret = ret ^ hash<Contact>()(*it);
				}
				return ret;
			}
		};
	}
}

#endif /* PSEUDOCONTACT_H_ */
