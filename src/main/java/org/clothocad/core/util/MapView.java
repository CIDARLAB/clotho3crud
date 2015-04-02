package org.clothocad.core.util;

import java.util.AbstractMap;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/** A read-only view of Map objects
 * Passes all operations through to the underlying source map, except for
 * mutator methods. Some methods such as entrySet() return a Set view, whose
 * mutators are also blocked. Mutators of the keys and values themselves are
 * not intercepted. Aside from the blocked mutators, MapView objects should
 * behave identically to their underlying source objects.
 */
public class MapView<K, V> implements Map<K, V> {
    private final Map<K, V> source;
    private MapView(final Map<K, V> source) {this.source = source;}

    @Override public void clear()
    { throw new UnsupportedOperationException(); }

    @Override public V put(final K x, final V y)
    { throw new UnsupportedOperationException(); }

    @Override public void putAll(final Map<? extends K, ? extends V> x)
    { throw new UnsupportedOperationException(); }

    @Override public V remove(final Object x)
    { throw new UnsupportedOperationException(); }

    @Override public boolean containsKey(final Object x)
    { return source.containsKey(x); }

    @Override public boolean containsValue(final Object x)
    { return source.containsValue(x); }

    @Override public Collection<V> values()
    { return new SetView<V>(source.values()); }

    @Override public Set<Map.Entry<K, V>> entrySet()
    { return new EntrySetView<K, V>(source.entrySet()); }

    @Override public boolean equals(final Object x) {return source.equals(x);}
    @Override public V get(final Object x)          {return source.get(x);}
    @Override public int hashCode()                 {return source.hashCode();}
    @Override public boolean isEmpty()              {return source.isEmpty();}
    @Override public int size()                     {return source.size();}
    @Override public Set<K> keySet() {return new SetView<K>(source.keySet());}

    /** Carefully wraps source in a MapView.
     * If source is null or is a MapView instance, returns source.
     * Otherwise, wraps source in a new MapView object.
     */
    public static <K, V> MapView<K, V> wrap(final Map<K, V> source) {
        if (source == null)
            return null;
        try {
            return (MapView<K, V>) source;
        } catch (ClassCastException e) {
            return new MapView<K, V>(source);
        }
    }

    private static abstract class SimpleIterator<E> implements Iterator<E> {
        @Override public void remove()
        { throw new UnsupportedOperationException(); }
    }

    private static class IteratorView<E> extends SimpleIterator<E> {
        private final Iterator<E> source;
        private IteratorView(final Iterator<E> source) {this.source = source;}
        @Override public boolean hasNext() {return source.hasNext();}
        @Override public E next() {return source.next();}
    }

    private static class EntrySetIteratorView<K, V>
    extends SimpleIterator<Map.Entry<K, V>> {
        private final Iterator<Map.Entry<K, V>> source;
        @Override public boolean hasNext() {return source.hasNext();}

        private EntrySetIteratorView(final Iterator<Map.Entry<K, V>> source)
        { this.source = source; }

        @Override public Map.Entry<K, V> next()
        { return new AbstractMap.SimpleImmutableEntry<K, V>(source.next()); }
    }

    private static abstract class AbstractSetView<E> implements Set<E> {
        abstract protected Collection<E> getSource();
        @Override public int hashCode() {return getSource().hashCode();}
        @Override public boolean isEmpty() {return getSource().isEmpty();}
        @Override public int size() {return getSource().size();}
        @Override public Object[] toArray() {return getSource().toArray();}

        @Override public <T> T[] toArray(final T[] x)
        { return getSource().toArray(x); }

        @Override public boolean add(final E x)
        { throw new UnsupportedOperationException(); }

        @Override public boolean addAll(final Collection<? extends E> x)
        { throw new UnsupportedOperationException(); }

        @Override public void clear()
        { throw new UnsupportedOperationException(); }

        @Override public boolean remove(final Object x)
        { throw new UnsupportedOperationException(); }

        @Override public boolean removeAll(final Collection<?> x)
        { throw new UnsupportedOperationException(); }

        @Override public boolean retainAll(final Collection<?> x)
        { throw new UnsupportedOperationException(); }

        @Override public boolean contains(final Object x)
        { return getSource().contains(x); }

        @Override public boolean containsAll(final Collection<?> x)
        { return getSource().containsAll(x); }

        @Override public boolean equals(final Object x)
        { return getSource().equals(x); }
    }

    private static class SetView<E> extends AbstractSetView<E> {
        private final Collection<E> source;
        private SetView(final Collection<E> source) {this.source = source;}
        @Override protected Collection<E> getSource() {return this.source;}

        @Override public Iterator<E> iterator()
        { return new IteratorView<E>(this.source.iterator()); }
    }

    private static class EntrySetView<K, V>
    extends AbstractSetView<Map.Entry<K, V>> {
        private final Set<Map.Entry<K, V>> source;
        @Override protected Collection<Map.Entry<K, V>> getSource()
        { return this.source; }

        private EntrySetView(final Set<Map.Entry<K, V>> source)
        { this.source = source; }

        @Override public Iterator<Map.Entry<K, V>> iterator()
        { return new EntrySetIteratorView<K, V>(this.source.iterator()); }
    }
}
