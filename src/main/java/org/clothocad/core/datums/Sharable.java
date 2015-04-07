package org.clothocad.core.datums;

import org.clothocad.model.Person;

public interface Sharable  {
    //Metadata for all Sharables
    public ObjectId getId();
    public Person getAuthor();
    public String getIcon();
    public String getName();
    public String getDescription();
}
