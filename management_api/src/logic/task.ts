import { uuid } from "../helpers/utils";
import { QUERY } from "../sql";
import { Task, NewTask } from "../types/task";

export async function getTasks(boardId: string): Promise<Task[]> {
	return QUERY.tasks.getTasks(boardId);
}

export async function createTask(task: NewTask): Promise<Task> {
	const newId = uuid();

	await QUERY.tasks.createTask({
		id: newId,
		status: task.done,
		boardId: task.boardId,
		title: task.title,
	});

	const createdTask = await QUERY.tasks.getTask(newId);
	if (!createdTask) throw "there is no task";

	return createdTask;
}

export async function updateTask(userId: string, taskId: string, newTask: {
	title: string,
	done: boolean,
	boardId: string,
}) {
	const task = await QUERY.tasks.getTask(taskId);
	if (!task) throw new Error('Task inexistente');

	const board = await QUERY.boards.getBoard(newTask.boardId);
	if (!board) throw new Error('Board inexistente');
	if (userId !== board?.ownerId) throw new Error("Sem acesso a essa board");
	if (board?.isArchived) throw new Error("A board está arquivada");

	await QUERY.tasks.updateTask(taskId, {
		title: newTask.title,
		done: newTask.done,
		boardId: newTask.boardId,
	});
}

export const tasks = {
	getTasks,
	createTask,
	updateTask,
};
