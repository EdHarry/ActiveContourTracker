/*
 Copied from Icy
 */
package activecontourtracker;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Processor class.<br>
 * Allow you to queue and execute tasks on a defined set of thread.
 *
 * @author stephane
 */
public
        class Processor extends ThreadPoolExecutor
{

    public static final
            int DEFAULT_MAX_WAITING = -1;
    public static final
            int DEFAULT_MAX_PROCESSING
                = Runtime.getRuntime().availableProcessors();

    protected
            class ProcessorThreadFactory implements ThreadFactory
    {

        String name;

        public
                ProcessorThreadFactory(String name)
        {
            super();

            SetName(name);
        }

        public
                String GetName()
        {
            return name;
        }

        public final
                void SetName(String value)
        {
            this.name = value;
        }

        String GetThreadName()
        {
            String result = name;

            return result;
        }

        @Override
        public
                Thread newThread(Runnable r)
        {
            final
                    Thread result = new Thread(r, GetThreadName());

            result.setPriority(priority);

            return result;
        }
    }

    protected
            class ProcessorRejectedExecutionHandler implements
            RejectedExecutionHandler
    {

        @Override
        public
                void rejectedExecution(Runnable r, ThreadPoolExecutor executor)
        {
            throw new RejectedExecutionException("Cannot add new task, "
                                                 + "ignore execution of " + r);
        }
    }

    protected
            class FutureTaskAdapter<T> extends FutureTask<T>
    {

        public
                Runnable runnable;
        public
                Callable<T> callable;
        final
                boolean handleException;

        public
                FutureTaskAdapter(Runnable runnable, T result, 
                                                     boolean handleException)
        {
            super(runnable, result);

            this.runnable = runnable;
            this.callable = null;
            this.handleException = handleException;
        }

        public
                FutureTaskAdapter(Runnable runnable, boolean handleException)
        {
            this(runnable, null, handleException);
        }

        public
                FutureTaskAdapter(Callable<T> callable, boolean handleException)
        {
            super(callable);

            this.runnable = null;
            this.callable = callable;
            this.handleException = handleException;
        }

        @Override
        protected
                void done()
        {
            super.done();

            if (handleException)
            {
                try
                {
                    get();
                }
                catch (InterruptedException | ExecutionException ex)
                {
                    Logger.getLogger(Processor.class.getName()).
                            log(Level.SEVERE, null, ex);
                }
            }
        }
    }

    /**
     * The minimum priority that a thread can have.
     */
    public final static
            int MIN_PRIORITY = Thread.MIN_PRIORITY;

    /**
     * The default priority that is assigned to a thread.
     */
    public final static
            int NORM_PRIORITY = Thread.NORM_PRIORITY;

    /**
     * The maximum priority that a thread can have.
     */
    public final static
            int MAX_PRIORITY = Thread.MAX_PRIORITY;

    /**
     * parameters
     */
    int priority;

    /**
     * internal
     */
    protected
            Runnable waitingExecution;
    protected
            long lastAdd;

    /**
     * Create a new Processor with specified number of maximum waiting and
     * processing tasks.<br>
     *
     * @param maxWaiting The length of waiting queue.
     * @param numThread The maximum number of processing thread.
     * @param priority Processor priority<br>
     * <code>Processor.MIN_PRIORITY</code><br>
     * <code>Processor.NORM_PRIORITY</code><br>
     * <code>Processor.MAX_PRIORITY</code>
     */
    public
            Processor(int maxWaiting, int numThread, int priority)
    {
        super(numThread, numThread, 2L, TimeUnit.SECONDS, (maxWaiting == -1)
                              ? new LinkedBlockingQueue<>() : 
                                       new LinkedBlockingQueue<>(maxWaiting));

        SetThreadFactory(new ProcessorThreadFactory("Processor"));
        SetRejectedExecutionHandler(new ProcessorRejectedExecutionHandler());
        AllowCoreThreadTimeOut(true);

        this.priority = priority;

        waitingExecution = null;
    }

    final
            void SetThreadFactory(ThreadFactory tf)
    {
        setThreadFactory(tf);
    }

    final
            void SetRejectedExecutionHandler(RejectedExecutionHandler handler)
    {
        setRejectedExecutionHandler(handler);
    }

    final
            void AllowCoreThreadTimeOut(boolean value)
    {
        allowCoreThreadTimeOut(value);
    }

    /**
     * Create a new Processor with specified number of maximum waiting and
     * processing tasks.
     *
     * @param maxWaiting The length of waiting queue.
     * @param numThread The maximum number of processing thread.
     */
    public
            Processor(int maxWaiting, int numThread)
    {
        this(maxWaiting, numThread, NORM_PRIORITY);
    }

    /**
     * Create a new Processor with specified number of processing thread.
     *
     * @param numThread The maximum number of processing thread.
     */
    public
            Processor(int numThread)
    {
        this(-1, numThread, NORM_PRIORITY);
    }

    /**
     * Create a new Processor with default number of maximum waiting and
     * processing tasks.
     */
    public
            Processor()
    {
        this(DEFAULT_MAX_WAITING, DEFAULT_MAX_PROCESSING);
    }

    @Override
    public
            boolean remove(Runnable task)
    {
        // don't forget to remove the reference here
        if (waitingExecution == task)
        {
            waitingExecution = null;
        }

        return super.remove(task);
    }

    /**
     * @param <T>
     * @param handledException if set to <code>true</code> then any occurring
     * exception during the runnable processing will be catch by
     * {@link IcyExceptionHandler}.
     * @param runnable the runnable task being wrapped
     * @param value the default value for the returned future
     * @return a <tt>RunnableFuture</tt> which when run will run the underlying
     * runnable and which, as a <tt>Future</tt>, will yield the given value as
     * its result and provide for cancellation of the underlying task.
     */
    protected
            <T> FutureTaskAdapter<T> NewTaskFor(boolean handledException,
                                                Runnable runnable, T value)
    {
        return new FutureTaskAdapter<>(runnable, value, handledException);
    }

    /**
     * @param <T>
     * @param handledException if set to <code>true</code> then any occurring
     * exception during the runnable processing will be catch by
     * {@link IcyExceptionHandler}.
     * @param callable the callable task being wrapped
     * @return a <tt>RunnableFuture</tt> which when run will call the underlying
     * callable and which, as a <tt>Future</tt>, will yield the callable's
     * result as its result and provide for cancellation of the underlying task.
     */
    protected
            <T> FutureTaskAdapter<T> NewTaskFor(boolean handledException,
                                                Callable<T> callable)
    {
        return new FutureTaskAdapter<>(callable, handledException);
    }

    @Override
    public
            void execute(Runnable task)
    {
        super.execute(task);
        // save the last executed task
        waitingExecution = task;
    }

    /**
     * Submit the given task (internal use only).
     *
     * @param <T>
     * @param task
     * @return
     */
    protected synchronized
            <T> FutureTask<T> Submit(FutureTaskAdapter<T> task)
    {
        execute(task);
        return task;
    }

    @Override
    public
            Future<?> submit(Runnable task)
    {
        if (task == null)
        {
            throw new NullPointerException();
        }

        return Submit(NewTaskFor(false, task, null));
    }

    @Override
    public
            <T> Future<T> submit(Runnable task, T result)
    {
        if (task == null)
        {
            throw new NullPointerException();
        }

        return Submit(NewTaskFor(false, task, result));
    }

    @Override
    public
            <T> Future<T> submit(Callable<T> task)
    {
        if (task == null)
        {
            throw new NullPointerException();
        }

        return Submit(NewTaskFor(false, task));
    }

    /**
     * Submits a Runnable task for execution and returns a Future representing
     * that task. The Future's <tt>get</tt> method will return <tt>null</tt>
     * upon <em>successful</em> completion.
     *
     * @param handleException if set to <code>true</code> then any occurring
     * exception during the runnable processing will be catch by
     * {@link IcyExceptionHandler}.
     * @param task the task to submit
     * @return a Future representing pending completion of the task
     * @throws RejectedExecutionException if the task cannot be scheduled for
     * execution
     * @throws NullPointerException if the task is null
     */
    public
            Future<?> Submit(boolean handleException, Runnable task)
    {
        if (task == null)
        {
            throw new NullPointerException();
        }

        return Submit(NewTaskFor(handleException, task, null));
    }

    /**
     * Submits a Runnable task for execution and returns a Future representing
     * that task. The Future's <tt>get</tt> method will return the given result
     * upon successful completion.
     *
     * @param <T>
     * @param handleException if set to <code>true</code> then any occurring
     * exception during the runnable processing will be catch by
     * {@link IcyExceptionHandler}.
     * @param task the task to submit
     * @param result the result to return
     * @return a Future representing pending completion of the task
     * @throws RejectedExecutionException if the task cannot be scheduled for
     * execution
     * @throws NullPointerException if the task is null
     */
    public
            <T> Future<T> Submit(boolean handleException,
                                 Runnable task, T result)
    {
        if (task == null)
        {
            throw new NullPointerException();
        }

        return Submit(NewTaskFor(handleException, task, result));
    }

    /**
     * Submits a value-returning task for execution and returns a Future
     * representing the pending results of the task. The Future's <tt>get</tt>
     * method will return the task's result upon successful completion.
     * <p>
     * If you would like to immediately block waiting for a task, you can use
     * constructions of the form
     * <tt>result = exec.submit(aCallable).get();</tt>
     * <p>
     * Note: The {@link Executors} class includes a set of methods that can
     * convert some other common closure-like objects, for example,
     * {@link java.security.PrivilegedAction} to {@link Callable} form so they
     * can be submitted.
     *
     * @param <T>
     * @param handleException if set to <code>true</code> then any occurring
     * exception during the runnable processing will be catch by
     * {@link IcyExceptionHandler}.
     * @param task the task to submit
     * @return a Future representing pending completion of the task
     * @throws RejectedExecutionException if the task cannot be scheduled for
     * execution
     * @throws NullPointerException if the task is null
     */
    public
            <T> Future<T> Submit(boolean handleException, Callable<T> task)
    {
        if (task == null)
        {
            throw new NullPointerException();
        }

        return Submit(NewTaskFor(handleException, task));
    }

    /**
     * Return true if one or more process are executing or we still have waiting
     * tasks.
     *
     * @return
     */
    public
            boolean IsProcessing()
    {
        return (getActiveCount() > 0) || HasWaitingTasks();
    }

    /**
     * Wait for all tasks completion
     */
    public
            void WaitAll()
    {
        while (IsProcessing())
        {
            Sleep(1);
        }
    }

    /**
     * shutdown and wait current tasks completion
     */
    public
            void ShutdownAndWait()
    {
        shutdown();
        while (!isTerminated())
        {
            Sleep(1);
        }
    }

    public static
            void Sleep(int milli)
    {
        try
        {
            Thread.sleep(milli);
        }
        catch (InterruptedException e)
        {
            // have to interrupt the thread
            Thread.currentThread().interrupt();
        }
    }

    /**
     * @return the priority
     */
    public
            int GetPriority()
    {
        return priority;
    }

    /**
     * @param priority the priority to set
     */
    public
            void SetPriority(int priority)
    {
        this.priority = priority;
    }

    /**
     * Return the thread name.
     *
     * @return
     */
    public
            String GetThreadName()
    {
        return ((ProcessorThreadFactory) getThreadFactory()).GetName();
    }

    /**
     * Set the wanted thread name.
     *
     * @param defaultThreadName
     */
    public
            void SetThreadName(String defaultThreadName)
    {
        ((ProcessorThreadFactory) getThreadFactory()).SetName(defaultThreadName);
    }

    /**
     * Get the number of free slot in queue
     *
     * @return
     */
    public
            int GetFreeSlotNumber()
    {
        return getQueue().remainingCapacity();
    }

    /**
     * Return true if queue is full
     *
     * @return
     */
    public
            boolean IsFull()
    {
        return GetFreeSlotNumber() == 0;
    }

    /**
     * Return waiting tasks
     *
     * @return
     */
    protected
            List<FutureTaskAdapter<?>> GetWaitingTasks()
    {
        final
                BlockingQueue<Runnable> q = getQueue();
        final
                List<FutureTaskAdapter<?>> result = new ArrayList<>();

        synchronized (q)
        {
            q.stream().forEach((r)
                    -> 
                    {
                        result.add((FutureTaskAdapter<?>) r);
            });
        }

        return result;
    }

    /**
     * Return waiting tasks for the specified Runnable instance
     *
     * @param task
     * @return
     */
    protected
            List<FutureTaskAdapter<?>> GetWaitingTasks(Runnable task)
    {
        final
                List<FutureTaskAdapter<?>> result = new ArrayList<>();

        // scan all tasks
        GetWaitingTasks().stream().filter((f)
                -> (f.runnable == task)).forEach((f)
                -> 
                {
                    result.add(f);
        });

        return result;
    }

    /**
     * Return waiting tasks for the specified Callable instance
     *
     * @param task
     * @return
     */
    protected
            List<FutureTaskAdapter<?>> GetWaitingTasks(Callable<?> task)
    {
        final
                List<FutureTaskAdapter<?>> result = new ArrayList<>();

        // scan all tasks
        GetWaitingTasks().stream().filter((f)
                -> (f.callable == task)).forEach((f)
                -> 
                {
                    result.add(f);
        });

        return result;
    }

    /**
     * Return the number of waiting task
     *
     * @return
     */
    public
            int GetWaitingTasksCount()
    {
        final
                int result = getQueue().size();

        // TODO : be sure that waitingExecution pass to null when task has been taken in account.
        // Queue can be empty right after a task submission.
        // For this particular case we return 1 if a task has been submitted
        // and not taken in account with a timeout of 1 second.
        if ((result == 0) && ((waitingExecution != null)
                              && ((System.currentTimeMillis() - lastAdd) < 1000)))
        {
            return 1;
        }

        return result;
    }

    /**
     * Return the number of task waiting in queue for the specified
     * <tt>Runnable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            int GetWaitingTasksCount(Runnable task)
    {
        int result = 0;

        result = GetWaitingTasks().stream().filter((f)
                -> (f.runnable == task)).map((_item) -> 1).
                reduce(result, Integer::sum);

        return result;
    }

    /**
     * Return the number of task waiting in queue for the specified
     * <tt>Callable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            int GetWaitingTasksCount(Callable<?> task)
    {
        int result = 0;

        result = GetWaitingTasks().stream().filter((f)
                -> (f.callable == task)).map((_item) -> 1).
                reduce(result, Integer::sum);

        return result;
    }

    /**
     * Return true if we have at least one task waiting in queue
     *
     * @return
     */
    public
            boolean HasWaitingTasks()
    {
        return (GetWaitingTasksCount() > 0);
    }

    /**
     * Return true if we have at least one task in queue for the specified
     * <tt>Runnable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            boolean HasWaitingTasks(Runnable task)
    {
        // scan all tasks

        return GetWaitingTasks().stream().anyMatch((f)
                -> (f.runnable == task));
    }

    /**
     * Return true if we have at least one task in queue for the specified
     * <tt>Callable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            boolean HasWaitingTasks(Callable<?> task)
    {
        // scan all tasks

        return GetWaitingTasks().stream().anyMatch((f)
                -> (f.callable == task));
    }

    /**
     * Remove first waiting task for the specified <tt>FutureTaskAdapter</tt>
     * instance.
     *
     * @param task
     * @return
     */
    protected
            boolean RemoveFirstWaitingTask(FutureTaskAdapter<?> task)
    {
        if (task == null)
        {
            return false;
        }

        synchronized (getQueue())
        {
            // remove first task of specified instance
            for (FutureTaskAdapter<?> f : GetWaitingTasks())
            {
                if (f == task)
                {
                    return remove(f);
                }
            }
        }

        return false;
    }

    /**
     * Remove first waiting task for the specified <tt>Runnable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            boolean ReoveFirstWaitingTask(Runnable task)
    {
        if (task == null)
        {
            return false;
        }

        synchronized (getQueue())
        {
            // remove first task of specified instance
            for (FutureTaskAdapter<?> f : GetWaitingTasks())
            {
                if (f.runnable == task)
                {
                    return remove(f);
                }
            }
        }

        return false;
    }

    /**
     * Remove first waiting task for the specified <tt>Callable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            boolean RemoveFirstWaitingTask(Callable<?> task)
    {
        if (task == null)
        {
            return false;
        }

        synchronized (getQueue())
        {
            // remove first task of specified instance
            for (FutureTaskAdapter<?> f : GetWaitingTasks())
            {
                if (f.callable == task)
                {
                    return remove(f);
                }
            }
        }

        return false;
    }

    /**
     * Remove all waiting tasks for the specified <tt>Runnable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            boolean RemoveWaitingTasks(Runnable task)
    {
        boolean result = false;

        synchronized (getQueue())
        {
            // remove all tasks of specified instance
            result = GetWaitingTasks(task).stream().map((f)
                    -> remove(f)).reduce(result, (accumulator, _item)
                                         -> accumulator | _item);
        }

        return result;
    }

    /**
     * Remove all waiting tasks for the specified <tt>Callable</tt> instance.
     *
     * @param task
     * @return
     */
    public
            boolean RemoveWaitingTasks(Callable<?> task)
    {
        boolean result = false;

        synchronized (getQueue())
        {
            // remove all tasks of specified instance
            result = GetWaitingTasks(task).stream().map((f)
                    -> remove(f)).reduce(result, (accumulator, _item)
                                         -> accumulator | _item);
        }

        return result;
    }

    /**
     * Clear all waiting tasks
     */
    public
            void RemoveAllWaitingTasks()
    {
        waitingExecution = null;

        synchronized (getQueue())
        {
            // remove all tasks
            getQueue().clear();
        }
    }

    @Override
    protected
            void beforeExecute(Thread t, Runnable r)
    {
        super.beforeExecute(t, r);

        // ok we can remove reference...
        waitingExecution = null;
    }
}
