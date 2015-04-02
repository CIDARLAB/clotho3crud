/*
Copyright (c) 2009 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */
package org.clothocad.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import lombok.Getter;
import lombok.Setter;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
 *
 * @author J. Christopher Anderson
 */
//name must be unique
public class Collection extends ObjBase {

    /**Constructor for collections from raw data
     *
     * @param collectionName = name of the Collection as a String
     * @param myauthor = author of Collection as a Person object
     */
    public Collection(String collectionName, String description, Person myauthor) {
        super(collectionName);
        this.description = description;
        this.author = myauthor;
        this.items = new HashMap<ObjBase,Object>();
    }

    /**Constructor for a transient Collection.  Transient collections
     * aren't saved to the database ever, but they can be passed around
     * via the Collector (they have UUIDs).
     *
     * @param collectionName = name of the Collection as a String
     * @param myauthor = author of Collection as a Person object
     */
    public Collection() {
        super("Results");
        description = "Transient Collection";
        items = new HashMap<ObjBase,Object>();
    }


    /**Abstract method addObject in general is for drag-and-drop events
     * Method gets called from the receiver of the drop
     *
     * @param dropObject is the object being dropped
     * @return is true if the drop is accepted for the receiver type
     */
    public boolean addObject(ObjBase dropObject, Object objitem) {
        
        if (items.containsKey(dropObject)) {
            return false;
        }

        AddAnyItem(dropObject,objitem);
        System.out.println("****Collection added an object " + dropObject.getName());
        return true;
    }

    /* SETTERS
     * */
    /**Remove an item from the Collection
     *
     * @param item = the item you want to remove
     */
    public boolean removeItem(ObjBase item) {
        if (!items.containsKey(item)){
            return false;
        }
        items.remove(item);
        return true;
    }

    /* GETTERS
     * */

    /**
     * Get all objects in the Collection
     * @return a List of ObjBase's in this Collection
     */
    public List<ObjBase> getAll() {
        List<ObjBase> out = new ArrayList<ObjBase>();
        for(Map.Entry<ObjBase, Object> entry: items.entrySet()){
            if (entry.getKey() == null) {
                continue;
            }
            out.add(entry.getKey());
        }
           
        /*for (ObjBase obj : items) {
            if (entry.getKey() == null) {
                continue;
            }
            out.add(obj);
        }*/
        return out;
    }

    /**
     * Get all ObjBase in this Collection and any ObjBases in a Collection
     * within this Collection
     * @param type the type of ObjBase you want
     * @return a HashSet full of ObjBases
     */
    public <T extends ObjBase> Set<T> recursiveGetAllOf(Class<T> type) {
        return recursiveRelay(type, new HashSet<ObjectId>());
    }

    <T extends ObjBase> Set<T> recursiveRelay(Class<T> type, Set<ObjectId> tested) {
        tested.add(getId());
        Set<T> out = new HashSet<T>();
        List<Collection> allmycollections = getAll(Collection.class);
        for (Collection childcoll : allmycollections) {
            if (tested.contains(childcoll.getId())) {
                continue;
            }
            Set<T> itscontents = childcoll.recursiveRelay(type, tested);
            for (T o : itscontents) {
                out.add(o);
            }
        }

        for (T o : getAll(type)) {
            out.add(o);
        }

        return out;
    }


    /**
     * Get a list of UUIDs for ObjBase in this Collection of a given type
     * @param type the ObjType desired
     * @return a Set of ObjBases
     */
    public <T extends ObjBase> Set<UUID> getAllLinksOf(Class<T> type) {
        HashSet out = new HashSet<UUID>();
        for(Map.Entry<ObjBase, Object> entry: items.entrySet())
        {
            if(type.isInstance(entry.getKey()))
            {
                out.add(entry.getKey().getId());
            }
        }
        return out;
    }

    /**Generic getter for all of whatever type
     *
     * @param type the ObjType desired
     * @return a List of ObjBase of that type
     */
    public <T extends ObjBase> List<T> getAll(Class<T> type) {
        List out = new ArrayList<T>();
        for(Map.Entry<ObjBase, Object> entry: items.entrySet())
        {
                if (type.isInstance(entry.getKey())) {
                out.add(entry.getKey());
            }
        }
        /*for (ObjBase item : items){
            if (type.isInstance(item)) {
                out.add(item);
            }
        }*/
        return out;
    }

    /**
     * This one probably shouldn't exists
     * deprecated
     * @param _myPart
     * @return
     *
    Deprecated
    public ArrayList<Plasmid> getPlasmidsOf(Part _myPart) {
        @SuppressWarnings(value = "unchecked")
        ArrayList<Plasmid> allplas = (ArrayList<Plasmid>) getAll(ObjType.PLASMID);
        ArrayList<Plasmid> out = new ArrayList<Plasmid>();
        for (Plasmid p : allplas) {
            if (p.getPart().getId().equals(_myPart.getId())) {
                out.add(p);
            }
        }
        return out;
    }*/

    /**
     * This one probably shouldn't exists
     * deprecated
     * @param _myPart
     * @return
     *
    Deprecated
    public ArrayList<PlasmidSample> getSamplesOf(Plasmid _myPlasmid) {
        @SuppressWarnings(value = "unchecked")
        ArrayList<Sample> allsam = (ArrayList<Sample>) getAll(ObjType.SAMPLE);
        ArrayList<PlasmidSample> out = new ArrayList<PlasmidSample>();
        for (Sample p : allsam) {
            PlasmidSample ps = (PlasmidSample) p;
            System.out.println("comparing " + ps.getPlasmid().getId() + "  " + _myPlasmid.getId());

            if (ps.getPlasmid().getId().equals(_myPlasmid.getId())) {
                out.add(ps);
            }
        }
        return out;
    }*/

    private void AddAnyItem(ObjBase item, Object obj) {
        if (!items.containsKey(item)){
            items.put(item, obj);
        }
    }

    public void removeAll() {
        items.clear();
    }

    public static Collection retrieveByName(String name) {
        throw new UnsupportedOperationException();
    }


    /*-----------------
    variables
    -----------------*/
    
    @Getter
    @Setter        
    private String description;
    
    @Getter
    @Reference
    private Person author;
    @ReferenceCollection
    //private List<ObjBase> items;
    private Map<ObjBase, Object> items;
            
   /* public static class CollectionDatum extends ObjBaseDatum {

        public Map<String, ObjType> uuidTypeHash = new HashMap<String, ObjType>();
        public Map<ObjType, HashSet<String>> typeUUIDHash = new EnumMap<ObjType, HashSet<String>>(ObjType.class);
        public ArrayList<String> itemUUIDs = new ArrayList<String>();
        public String _authorUUID;
        public String _description;

        @Override
        public ObjType getType() {
            return ObjType.COLLECTION;
        }
    }*/
}
