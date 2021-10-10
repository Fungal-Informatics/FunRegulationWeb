import { zap } from "../helpers/sql";

async function createUser(args: {
	id: string;
	name: string;
	email: string;
	password: string;
}) {
	const { id, name, email, password } = args;
	await zap.insert("users", {
		id,
		password,
		email,
		name,
		createdAt: new Date(),
		updatedAt: new Date(),
		isArchived: false,
	});
}

async function getUserByEmail(email: string) {
	return zap.selectOne("users", { email });
}

async function getUser(id: string) {
	return zap.selectOne("users", { id });
}

async function getUserByEmailAndPassword(email: string, password: string) {
	return zap.selectOne("users", { email, password });
}

async function setUserSessionKey(userId: string, key: string) {
	const currentSession = await zap.selectOne("sessions", { key });
	if (!currentSession) {
		await zap.insert("sessions", { userId, key });
	} else {
		await zap.update("sessions", { key }, { userId });
	}
}

async function removeUserSessionKey(key: string) {
	await zap.deletes("sessions", { key });
}

async function getUserIdByDeviceId(key: string) {
	const user = await zap.selectOne("sessions", { key });
	return user?.userId ?? null;
}

export const users = {
	createUser,
	getUserByEmail,
	getUser,
	getUserByEmailAndPassword,
	setUserSessionKey,
	getUserIdByDeviceId,
	removeUserSessionKey,
};
