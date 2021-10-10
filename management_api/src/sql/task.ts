import { zap } from "../helpers/sql";
import { Task } from "../types/task";

async function getTasks(boardId: string): Promise<Task[]> {
	return zap.select("tasks", { boardId });
}

async function createTask(task: {
	id: string;
	boardId: string;
	title: string;
	status: boolean;
}) {
	await zap.insert("tasks", {
		...task,
		done: task.status,
		createdAt: new Date(),
		isArchived: false,
		updatedAt: new Date(),
	});
}

async function getTask(taskId: string): Promise<Task | undefined> {
	return zap.selectOne("tasks", { id: taskId });
}

async function changeTaskStatus(
	taskId: string,
	status: boolean,
): Promise<void> {
	await zap.update(
		"tasks",
		{
			id: taskId,
		},
		{
			done: status,
			updatedAt: new Date(),
		},
	);
}

async function updateTask(taskId: string, newTask: {
	title: string,
	done: boolean,
	boardId: string,
}) {

	await zap.update(
		"tasks",
		{ ...newTask, updatedAt: new Date() },
		{ id: taskId },
	);
}

export const tasks = {
	getTasks,
	createTask,
	getTask,
	changeTaskStatus,
	updateTask,
};
