// Copied from Icy
package activecontourtracker;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

public
        class Tracking
{

    public static
            class Detection implements Cloneable
    {

        Vector2 position;
        int t;

        public
                Detection(Vector2 position, int t)
        {
            this.position = position;
            this.t = t;
        }

        @Override
        protected
                Detection clone() throws CloneNotSupportedException
        {
            super.clone();
            return new Detection(position.Copy(), t);
        }
    }

    public static
            class TrackGroup
    {

        private
                ArrayList<TrackSegment> trackSegmentList;

        public
                TrackGroup()
        {
            trackSegmentList = new ArrayList<>();
        }

        public
                ArrayList<TrackSegment> GetTrackSegmentList()
        {
            return trackSegmentList;
        }

        public
                void AddTrackSegment(TrackSegment ts)
        {
            if (ts.GetOwnerTrackGroup() != null)
            {
                System.err.println("The trackSegment is already owned "
                                   + "by another TrackGroup.");
                return;
            }

            ts.SetOwnerTrackGroup(this);
            trackSegmentList.add(ts);
        }

        public
                TrackSegment GetTrackSegmentWithDetection(Detection detection)
        {
            ArrayList<TrackSegment> trackSegList = GetTrackSegmentList();

            for (TrackSegment ts : trackSegList)
            {
                if (ts.ContainsDetection(detection))
                {
                    return ts;
                }
            }

            return null;
        }

        public
                void ClearAllTrackSegment()
        {
            ArrayList<TrackSegment> trackSegmentListCopy
                                    = new ArrayList<>(trackSegmentList);

            trackSegmentListCopy.stream().forEach((ts)
                    ->
            {
                RemoveTrackSegment(ts);
            });
        }

        public
                void RemoveTrackSegment(TrackSegment ts)
        {
            // should remove links
            ts.SetOwnerTrackGroup(null);
            ts.RemoveAllLinks();
            trackSegmentList.remove(ts);
        }
    }

    public static
            class TrackSegment implements Cloneable
    {

        private
                ArrayList<Detection> detectionList = new ArrayList<>();
        ArrayList<TrackSegment> previousList = new ArrayList<>();
        ArrayList<TrackSegment> nextList = new ArrayList<>();
        int id;

        public static
                HashMap<TrackSegment, Integer> idHashMapList = new HashMap<>(); // 1-1 hashmap.

        private final static
                AtomicInteger ID_GENERATOR = new AtomicInteger();

        public static
                TrackSegment GetTrackSegmentById(int id)
        {

            for (TrackSegment ts : idHashMapList.keySet())
            {
                if (idHashMapList.get(ts) == id)
                {
                    return ts;
                }
            }

            return null;
        }

        @Override
        protected
                Object clone() throws CloneNotSupportedException
        {

            TrackSegment cloneSegment = (TrackSegment) super.clone();

            for (Detection detection : detectionList)
            {
                Detection cloneDetection = (Detection) detection.clone();
                //cloneDetection.setOwnerTrackSegment(cloneSegment);
                cloneSegment.detectionList.add(cloneDetection);
            }

            previousList = new ArrayList<>(previousList);
            nextList = new ArrayList<>(nextList);

            GenerateId();

            return cloneSegment;

        }

        public
                TrackSegment()
        {
            GenerateId();
        }

        private

                void GenerateId()
        {

            int localId = -1;
            
            idHashMapList.remove(this);
            boolean passed = false;
            while (!passed)
            {
                synchronized (ID_GENERATOR)
                {
                    //this.id = (int) (Math.random() * Integer.MAX_VALUE);
                    
                    localId = ID_GENERATOR.getAndIncrement();
                    
                    if (idHashMapList.containsValue(localId) == false)
                    {
                        passed = true;
                    }
                }
            }
            
            SetId(localId);
        }

        /**
         * Constructor with a list of detection
         */
        TrackSegment(ArrayList<Detection> detectionList)
        {
            this.detectionList = detectionList;
            GenerateId();
        }

        public
                int GetId()
        {
            return id;
        }

        public
                void SetId(int id)
        {
            if (idHashMapList.containsValue(id) == false)
            {
                this.id = id;
                idHashMapList.put(this, id);
            }
            else
            {
                System.out.println("track id already loaded");
            }
        }

        public
                boolean ContainsDetection(Detection detection)
        {
            return detectionList.contains(detection);
        }

        /**
         * Add a detection in the segmentTrack
         */
        void AddDetection(Detection detection)
        {
            if (detectionList.size() > 0)
            {
                Detection detectionPrevious = GetLastDetection();
                if (detection.t != detectionPrevious.t + 1)
                {
                    System.err.println("TrackSegment : The detection must be "
                                       + "added with consecutive T value. Detection "
                                       + "was not added");
                    return;
                }
            }
            detectionList.add(detection);
        }

        public
                void RemoveDetection(Detection detection)
        {
            detectionList.remove(detection);
        }

        /**
         * Remove a detection in the segmentTrack
         */
        void RemoveLastDetection()
        {
            detectionList.remove(GetLastDetection());
        }

        /**
         * Add a TrackSegment before this trackSegment
         */
        void AddPrevious(TrackSegment trackSegment)
        {
            previousList.add(trackSegment);
            trackSegment.nextList.add(this);
        }

        /**
         * Remove a TrackSegment before this trackSegment
         */
        void RemovePrevious(TrackSegment trackSegment)
        {
            previousList.remove(trackSegment);
            trackSegment.nextList.remove(this);
        }

        /**
         * Add a TrackSegment after this trackSegment
         */
        void AddNext(TrackSegment trackSegment)
        {
            nextList.add(trackSegment);
            trackSegment.previousList.add(this);
        }

        /**
         * Remove a TrackSegment after this trackSegment
         */
        void RemoveNext(TrackSegment trackSegment)
        {
            nextList.remove(trackSegment);
            trackSegment.previousList.remove(this);
        }

        /**
         * return first detection ( should be first in time too )
         */
        Detection GetFirstDetection()
        {
            if (detectionList.isEmpty())
            {
                return null;
            }

            return detectionList.get(0);
        }

        /**
         * return detection at index i
         */
        Detection GetDetectionAt(int i)
        {
            return detectionList.get(i);
        }

        /**
         * return detection at time t
         */
        Detection GetDetectionAtTime(int t)
        {
            for (Detection detection : detectionList)
            {
                if (detection.t == t)
                {
                    return detection;
                }
            }
            return null;
        }

        /**
         * return detection list WARNING: User should use addDetection and
         * removeDetection instead of doing it himself using direct access to
         * the ArrayList. Using addDetection or removeDetection ensure the
         * property ownerTrackSegment of the Detection to be correctly updated.
         */
        ArrayList<Detection> GetDetectionList()
        {
            return detectionList;
        }

        /**
         * return last detection ( should be last in time too )
         */
        Detection GetLastDetection()
        {
            if (detectionList.isEmpty())
            {
                return null;
            }

            return detectionList.get(detectionList.size() - 1);
        }

        /**
         * return detection index
         */
        int GetDetectionIndex(Detection detection)
        {
            return detectionList.indexOf(detection);
        }

        TrackGroup ownerTrackGroup = null;

        public
                void SetOwnerTrackGroup(TrackGroup tg)
        {

            ownerTrackGroup = tg;

        }

        public
                TrackGroup GetOwnerTrackGroup()
        {
            return ownerTrackGroup;
        }

        public
                void RemoveId()
        {

            idHashMapList.remove(this);

        }

        public
                void RemoveAllLinks()
        {

            ArrayList<TrackSegment> previousListCopy
                                    = new ArrayList<>(previousList);
            previousListCopy.stream().forEach(this::RemovePrevious);

            ArrayList<TrackSegment> nextListCopy = new ArrayList<>(nextList);
            nextListCopy.stream().forEach(this::RemoveNext);

        }

    }
}
