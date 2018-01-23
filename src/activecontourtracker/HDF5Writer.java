/*

Copied from openanalytics/phaedra BaseHDF5File.java
https://github.com/openanalytics/phaedra/blob/master/eu.openanalytics.phaedra.base.hdf5/src/eu/openanalytics/phaedra/base/hdf5/BaseHDF5File.java

 */
package activecontourtracker;

import java.io.IOException;
import ch.systemsx.cisd.hdf5.HDF5CompoundMemberInformation;
import ch.systemsx.cisd.hdf5.HDF5CompoundType;
import ch.systemsx.cisd.hdf5.HDF5DataClass;
import ch.systemsx.cisd.hdf5.HDF5DataSetInformation;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ncsa.hdf.hdf5lib.exceptions.HDF5JavaException;

class HDF5Writer
{

    private final
            String path;
    private
            IHDF5Writer writer;

    public
            HDF5Writer(String path)
    {
        this.path = path;
        Open();
    }

    private
            void Open()
    {
        if (writer == null)
        {
            this.writer = HDF5Factory.open(path);
        }
    }

    void Close()
    {
        if (writer != null)
        {
            try
            {
                writer.flushSyncBlocking();
            }
            catch (Throwable t)
            {
                // Ignore errors if file is already closed.
            }

            writer.close();
        }

        writer = null;
    }

    private
            long[] GetDataDimensions(String path)
    {
        HDF5DataSetInformation info = writer.getDataSetInformation(ResolvePath(path));
        if (info == null)
        {
            return new long[0];
        }
        return info.getDimensions();
    }

    void WriteCompoundData(String path, String[] columnNames, Object[][] data) throws IOException
    {
        WriteCompoundData(path, columnNames, data, false);
    }

    void WriteCompoundData(String path, String[] columnNames, Object[][] data, boolean overwrite) throws IOException
    {
        int columns = columnNames.length;
        Object[] dataTypes = new Object[columns];

        for (int col = 0;
             col < columns;
             ++col)
        {
            Object colSample = null;
            // Keep looping the rows in a column till a value is found.
            int row = 0;
            for (;
                row < data.length;
                ++row)
            {
                colSample = data[row][col];
                if (colSample != null)
                {
                    break;
                }
            }
            if (colSample == null)
            {
                throw new IOException("Unsupported data class 'null' for column " + columnNames[col] + " in data " + path);
            }
            if (colSample instanceof Number)
            {
                if (colSample instanceof Long)
                {
                    dataTypes[col] = Long.MAX_VALUE;
                }
                else if (colSample instanceof Integer)
                {
                    dataTypes[col] = Integer.MAX_VALUE;
                }
                 else if (colSample instanceof Double)
                {
                    dataTypes[col] = Double.NaN;
                }
                else
                {
                    dataTypes[col] = Float.NaN;
                }
            }
            if (colSample instanceof String)
            {
                String longest = colSample.toString();
                for (;
                    row < data.length;
                    ++row)
                {
                    colSample = data[row][col];
                    if (colSample != null)
                    {
                        String newString = colSample.toString();
                        if (newString.length() > longest.length())
                        {
                            longest = newString;
                        }
                    }
                }
                dataTypes[col] = longest;
            }
        }

        WriteCompoundData(path, columnNames, dataTypes, data, overwrite);
    }

    void WriteCompoundData(String path, String[] columnNames, Object[] dataTypes, Object[][] data) throws IOException
    {
        WriteCompoundData(path, columnNames, dataTypes, data, false);
    }

    void WriteCompoundData(String path, String[] columnNames, Object[] dataTypes, Object[][] data, boolean overwrite) throws IOException
    {
        try
        {
            // Create the column data types.
            HDF5CompoundType<Object[]> inferredType = writer.compounds().getInferredType(columnNames, dataTypes);
            // Check if the Compound Table already exists.
            if (!writer.exists(path))
            {
                // Compound Table does not exist yet. Create it.
                writer.compounds().createArray(path, inferredType, 0, HDF5GenericStorageFeatures.GENERIC_CHUNKED_KEEP);
                writer.compounds().writeArray(path, inferredType, data, HDF5GenericStorageFeatures.GENERIC_CHUNKED_KEEP);
            }
            else
            {
                // Compound Table does exist. See if the same columns are present.
                HDF5CompoundMemberInformation[] dataSetInfos = writer.compounds().getDataSetInfo(path);
                HDF5CompoundMemberInformation[] newDataSetInfos = inferredType.getCompoundMemberInformation();
                if (overwrite || !IsIdenticalDataTypes(dataSetInfos, newDataSetInfos))
                {
                    // Delete previous compound & overwrite existing data.
                    writer.delete(path);
                    writer.compounds().createArray(path, inferredType, 0, HDF5GenericStorageFeatures.GENERIC_CHUNKED_KEEP);
                    writer.compounds().writeArray(path, inferredType, data, HDF5GenericStorageFeatures.GENERIC_CHUNKED_KEEP);
                    // TODO: Previous Compound Types remain in the HDF5 file. Find a proper way to delete them.
                }
                else
                {
                    // Get the size for the offset.
                    long[] dims = GetDataDimensions(path);
                    writer.compounds().writeArrayBlockWithOffset(path, inferredType, data, dims[0]);
                }
            }
            writer.flush();
        }
        catch (HDF5JavaException e)
        {
            throw new IOException("Failed to write data " + path, e);
        }
    }

    private
            String ResolvePath(String path)
    {
        if (writer.exists(path))
        {
            // Resolve references, if any.
            if (writer.isDataSet(path))
            {
                HDF5DataClass dataClass = writer.getDataSetInformation(path).getTypeInformation().getDataClass();
                if (dataClass.toString().equals("REFERENCE"))
                {
                    return writer.readObjectReference(path);
                }
            }
            return path;
        }

        return null;
    }

    private
            boolean IsIdenticalDataTypes(HDF5CompoundMemberInformation[] dataSetInfos, HDF5CompoundMemberInformation[] newDataSetInfos)
    {
        boolean isIdentical = dataSetInfos.length == newDataSetInfos.length;
        if (isIdentical)
        {
            for (int i = 0;
                 i < dataSetInfos.length;
                 ++i)
            {
                if (!dataSetInfos[i].equals(newDataSetInfos[i]))
                {
                    isIdentical = false;
                    break;
                }
            }
        }
        return isIdentical;
    }

}
