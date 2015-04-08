package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
 *
 * @author J. Christopher Anderson
 */
@NoArgsConstructor
public class Person extends SharableObjBase {

    @Getter
    @Setter
    @Reference
    private Collection myCollection;

    @Getter
    @Setter
    private String givenName, surName, nickName, emailAddress, snailMailAddress;

    @Getter
    @Setter
    private String primaryEmail;

    //@Getter
    //@Setter
    //private boolean isPrimaryAccount;

    /**Constructor from raw data
     *
     * @param displayname = String of author such as "JCAnderson"
     * @param affiliation = String of affiliation such as "UC Berkeley"
     */
    //unique name criterion
    //valid or nonexistent email
    public Person( String displayname) {
        //XXX:  Do people have authors?
        super(displayname, null);
        //changePassword( rawPassword );
        myCollection = new Collection();
        //biography = new WikiText("");
    }

    /* SETTERS
     * */

    /**
     * Public accessible method for setting a Person as administrator.
     *
     * An administrator must be the current user for this to do anything.
     * -OR-
     * If the database has only 3 people in it the change is also allowed
     *
     * Admins have the ability to clear a person's password and possibly other
     * sensitive things.
     *
     * @param isit
     */
    public final void setAsAdministrator(boolean isit) {
        throw new UnsupportedOperationException("Not implemented yet.");
    }

    /* SETTERS
     * */

    /**
     * Method for clearing the password so that a new one can be entered
     * An administrator must be logged in to use this
     */
    /*
    public final boolean clearPassword() {
        throw new UnsupportedOperationException("Not implemented yet.");
    }*/

    /**
     * Plugin-accessible method for changing the Person's password
     * @param raw
     */
    /*
    public final void changePassword( String raw ) {
        throw new UnsupportedOperationException("Not implemented yet.");
    }*/

    /**
     * Check the Person's password
     * @return true if the user successfully provided the correct password
     */
    /*
    public final boolean checkPassword() {
        throw new UnsupportedOperationException("Not implemented yet.");
    }*/

    /**
     * Login this Person.  This involves validating that the password has been confirmed.
     */
    /*
    public final void login() {
        throw new UnsupportedOperationException("Not implemented yet.");
    }*/

    /**
     * Plugin-accessible call to determine if the Person is logged in
     * @return true if this Person is logged in
     */
    
    /*
    public final boolean isLoggedIn() {
        throw new UnsupportedOperationException();
    }*/

    /**
     * Change the User's first name
     * @param str a String
     */
    /*public final void changeGivenName( String str ) {
        if(!hasChangeClearance()) {
            //fireData(new RefreshEvent(this, RefreshEvent.Condition.GIVEN_NAME_CHANGED));
            return;
        }
        //addUndo( "_given_name", _personDatum._given_name, str );
        _personDatum._given_name = str;
        //setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.GIVEN_NAME_CHANGED);
    }*/

    /**
     * Change the user's last name (surname)
     * @param str a String
     */
   /* public final void changeSurName( String str ) {
        if(!hasChangeClearance()) {
            fireData(new RefreshEvent(this, RefreshEvent.Condition.SURNAME_CHANGED));
            return;
        }
        addUndo( "_surname", _personDatum._surname, str );
        _personDatum._surname = str;
        setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.SURNAME_CHANGED);
    }

    /**
     * Change the name a user wishes to be called
     * @param str a String
     */
    /*public final void changeNickName( String str ) {
        if(!hasChangeClearance()) {
            fireData(new RefreshEvent(this, RefreshEvent.Condition.NICKNAME_CHANGED));
            return;
        }
        addUndo( "_nick_name", _personDatum._nick_name, str );
        _personDatum._nick_name = str;
        setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.NICKNAME_CHANGED);
    }



    /**
     * Change the user's display name
     * @param str a String
     */
    public final void setDisplayName(String str) {
        setName(str);
    }

    /* GETTERS
     * */

    /**
     * Retrieve a Person ObjBase from the database using their display name.
     * This will only query Persons saved to the database, not local ones.
     * @param name
     * @return
     */
    public static Person retrieveByName(String name) {
        throw new UnsupportedOperationException();
    }

    //Removing this since Person should not have password
    /**
     * It checks a password to see if it matches the user's password and returns
     * true if they match, otherwise returns false.
     * It only stores the encrypted version of the password using SHA-1 hashing.
     *
     * @param raw the raw password supplied by user
     * @return true if it's a match
     */
    /*public final boolean checkPassword( String raw ) {
     throw new UnsupportedOperationException();
    }*/

    /*private boolean hasChangeClearance() {
        //If it's a new Person, changes are OK
        if(_isBrandNew) {
            return true;
        }
        Person user = Collector.getCurrentUser();
        //If nobody's logged in, they aren't allowed to change things
        if(user==null) {
            JOptionPane.showMessageDialog( null, "You aren't logged in.  You aren't allowed to change this data.", "Forbidden data change", JOptionPane.OK_OPTION );
            return false;
        }
        //If a person is editing their own fields, that's ok
        if(user.getUUID().equals(this.getUUID())) {
            return true;
        }
        //If the logged in user is an admin, they can edit fields
        if(user._personDatum._isAdministrator) {
            return true;
        }
        JOptionPane.showMessageDialog( null, "Only the user herself or admins can change this data.", "Forbidden data change", JOptionPane.OK_OPTION );
        return false;
    }*/

    /**
     * Is the person an administrator?
     * @return true if they are
     */
   // public boolean isAdmin() {
   //     return _personDatum._isAdministrator;
   // }

    /**
     * Get the person's login name
     * @return a String
     */
    public String getDisplayName() {
        return getName();
    }

    /**
     * Get the personal Collection of this object
     * @return a Collection ObjBase
     */
    public Collection getHerCollection() {
        return myCollection;
    }

    @Override
    public String toString() {
        return getName();
    }

}
