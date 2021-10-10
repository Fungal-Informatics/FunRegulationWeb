import { LOGIC } from "../logic";
import { api } from "../sdkgen/api-generated";

api.fn.getTasks = async (ctx, { boardId }) => {
    return LOGIC.tasks.getTasks(boardId);
};

api.fn.createTask = async (ctx, { boardId, task }) => {
    return LOGIC.tasks.createTask({
        boardId,
        ...task,
    })
};

api.fn.updateTask = async (ctx, {
    taskId,
    task,
}) => {
    await LOGIC.tasks.updateTask(ctx.request.extra.userId, taskId, task);
};

api.fn.archiveTask;
