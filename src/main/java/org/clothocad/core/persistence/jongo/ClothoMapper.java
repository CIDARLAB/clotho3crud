/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.ObjectMapper;
import static com.fasterxml.jackson.databind.SerializationFeature.FAIL_ON_EMPTY_BEANS;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.jackson.IdRenamingMixin;
import org.clothocad.core.persistence.jongo.ClothoObjectIdUpdater.ClothoIdFieldSelector;
import org.clothocad.core.util.JSON.ClothoJacksonModule;
import org.jongo.Mapper;
import org.jongo.ObjectIdUpdater;
import org.jongo.marshall.Marshaller;
import org.jongo.marshall.jackson.configuration.AbstractMappingBuilder;
import org.jongo.marshall.jackson.configuration.MapperModifier;
import org.jongo.marshall.jackson.configuration.Mapping;
import org.jongo.query.BsonQueryFactory;
import org.jongo.query.QueryFactory;

/**
 *
 * @author spaige
 */
public class ClothoMapper implements Mapper {
    
    private final RefResolvingJacksonEngine engine;
    private final ObjectIdUpdater objectIdUpdater;
    private final QueryFactory queryFactory;

    public ClothoMapper(){
        //behold my mistreatment of the builder pattern and weep
        ClothoMapperBuilder builder = new ClothoMapperBuilder();
        Mapping mapping = builder.createMapping();
        engine = new RefResolvingJacksonEngine(mapping);
        objectIdUpdater = new ClothoObjectIdUpdater(new ClothoIdFieldSelector());
        queryFactory = new BsonQueryFactory(engine);
    }

    @Override
    public Marshaller getMarshaller() {
        return engine;
    }

    @Override
    public ExtendedUnmarshaller getUnmarshaller() {
        return engine;
    }

    @Override
    public ObjectIdUpdater getObjectIdUpdater() {
        return objectIdUpdater;
    }

    @Override
    public QueryFactory getQueryFactory() {
        return queryFactory;
    }
    
    //exposing some protected utility methods
    public class ClothoMapperBuilder extends AbstractMappingBuilder<ClothoMapperBuilder>{

        public ClothoMapperBuilder() {
            super();
        }
        
        @Override
        public Mapping createMapping(){
            //add customizations
            addModifier(new MapperModifier() {
                @Override
                public void modify(ObjectMapper mapper) {
                        mapper.disable(FAIL_ON_EMPTY_BEANS);
                        //write types into serialized objects
                        //mapper.enableDefaultTyping(ObjectMapper.DefaultTyping.NON_CONCRETE_AND_ARRAYS);
                        mapper.disableDefaultTyping();
                        //mapper.setDefaultTyping(typer);
                        mapper.registerModule(new ClothoJacksonModule());
                        //add "id" -> "_id" renaming mixin
                        mapper.addMixInAnnotations(ObjBase.class, IdRenamingMixin.class);
                    }
                });
            return super.createMapping();
        }

        @Override
        protected ClothoMapperBuilder getBuilderInstance() {
            return this;
        }
    }
}
