package org.clothocad.core.communication;

import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.Map;
import org.clothocad.core.util.MapView;

/**
 * @author spaige
 */
public class Message {
    private final Channel channel;
    private final Object data;
    private final String requestId;
    private final Map<MessageOption, Object> options;

    @JsonProperty("channel")
    public Channel getChannel() { return channel; }

    @JsonProperty("data")
    public Object getData() { return data; }

    @JsonProperty("requestId")
    public String getRequestId() { return requestId; }

    @JsonProperty("options")
    public Map<MessageOption, Object> getOptions() { return options; }

    public Message(Channel channel, Object data, String requestId){
        this(channel, data, requestId, null);
    }
    
    public Message(@JsonProperty("channel") Channel channel,
                   @JsonProperty("data") Object data,
                   @JsonProperty("requestId") String requestId,
                   @JsonProperty("options") Map<MessageOption, Object> options) {
        this.channel = channel;
        this.data = data;
        this.requestId = requestId;
        this.options = MapView.wrap(options);
    }
}
