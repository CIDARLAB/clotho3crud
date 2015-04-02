package org.clothocad.core.communication;

public enum Channel {

//HUMAN INTERACTION
    autocomplete, //Return potential commands that start with this substring
    autocompleteDetail, //Return the metadata for a Sharable
    recent,  //Return N most recently used sharables
    submit, //Run this sloppy or concrete command
    clear, //Replace the scriptengine behind mind with a fresh version
    createUser, //Create a new user account for me.
    linkPerson, //Link an existing? or new Person to an existing Person Object
    getAssociatedPerson, //Get a list of all Person Objects associated with that primaryEmail.
    learn, //associate this String with this execution statement
    like,
    dislike,

    

//LOGGING AND MESSAGING
    log, //write a message to the console.log spot
    say, //Post a message to all search bars for this Person
    note, //Remember this message in my personal log (not implemented)
    alert, //Post an alert to this Person, wherever they say they want to be contacted in preferences.  This isn't a javascript alert().  That is just one implementation.

    
//DATA MANIPULATION
    get, //Get the Sharable with this uuid
    set, //Change the state of the Sharable with this uuid with this new JSON data
    create, //Create a new sharable with this Schema
    destroy, //Delete a sharable with this uuid
    query, //Find all Sharables that satisfy this formally-expressed constraint
    
    convert, //Convert one data object from one type to another. 
    
    getAll,
    createAll,
    destroyAll,
    setAll,
    queryOne,
    
    validate,
    
//EXECUTION
    run, //run this Function
    listen, //Listen for events, and in response do this execution statement
    unlisten, //Remove a listener for an event
    
//SECURITY
    grant,
    grantAll,
    login, //Log me into Clotho on this client with this login/password
    logout, //Log me out of Clotho
    changePassword, //Change my password to this new value
    createAccount,
    updatePassword, //Update the userPassword
    
//OTHER
    reloadModels, //convenience function for reloading test data
}
